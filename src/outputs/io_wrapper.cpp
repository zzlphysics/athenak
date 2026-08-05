//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file io_wrapper.cpp
//! \brief functions that provide wrapper for MPI-IO versus serial input/output

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "io_wrapper.hpp"

namespace {

[[noreturn]] void IOWrapperFatal(const std::string &message) {
  std::cerr << "### FATAL I/O ERROR";
#if MPI_PARALLEL_ENABLED
  std::cerr << " on rank " << global_variable::my_rank;
#endif
  std::cerr << ": " << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

std::size_t DataSize(const std::string &datatype) {
  if (datatype == "byte") return sizeof(char);
  if (datatype == "int") return sizeof(int);
  if (datatype == "float") return sizeof(float);
  if (datatype == "double") return sizeof(double);
  if (datatype == "Real") return sizeof(Real);
  IOWrapperFatal("unrecognized datatype '" + datatype + "'");
}

void SeekStdioForWrite(FILE *file, IOWrapperSizeT offset) {
  if (std::fseek(file, offset, SEEK_SET) != 0) {
    IOWrapperFatal("could not seek output file: " + std::string(std::strerror(errno)));
  }
}

std::size_t WriteStdio(const void *buf, IOWrapperSizeT count,
                       const std::string &datatype, FILE *file) {
  const std::size_t written = std::fwrite(buf, DataSize(datatype), count, file);
  if (written != count) {
    std::ostringstream msg;
    msg << "short write: expected " << count << " elements, wrote " << written;
    if (errno != 0) msg << " (" << std::strerror(errno) << ")";
    IOWrapperFatal(msg.str());
  }
  return written;
}

#if MPI_PARALLEL_ENABLED
std::string MPIErrorMessage(const char *operation, int error_code) {
  char mpi_message[MPI_MAX_ERROR_STRING];
  int message_length = 0;
  MPI_Error_string(error_code, mpi_message, &message_length);
  return std::string(operation) + ": " +
         std::string(mpi_message, static_cast<std::size_t>(message_length));
}

MPI_Datatype MPIDataType(const std::string &datatype) {
  if (datatype == "byte") return MPI_BYTE;
  if (datatype == "int") return MPI_INT;
  if (datatype == "float") return MPI_FLOAT;
  if (datatype == "double") return MPI_DOUBLE;
  if (datatype == "Real") return MPI_ATHENA_REAL;
  IOWrapperFatal("unrecognized datatype '" + datatype + "'");
}

std::size_t CheckMPIWrite(const MPI_Status &status, MPI_Datatype datatype,
                          IOWrapperSizeT expected, const char *operation) {
  int written = MPI_UNDEFINED;
  const int error_code = MPI_Get_count(&status, datatype, &written);
  if (error_code != MPI_SUCCESS) {
    IOWrapperFatal(MPIErrorMessage("MPI_Get_count after output", error_code));
  }
  if (written == MPI_UNDEFINED || static_cast<IOWrapperSizeT>(written) != expected) {
    std::ostringstream msg;
    msg << operation << " reported a short or undefined write: expected " << expected
        << " elements, wrote " << written;
    IOWrapperFatal(msg.str());
  }
  return static_cast<std::size_t>(written);
}

std::size_t ReadPerRankFile(void *buf, std::size_t size, std::size_t count,
                            FILE *file) {
  const std::size_t nread = std::fread(buf, size, count, file);
  if (nread != count) {
    std::cerr << "Per-rank file read was truncated on rank "
              << global_variable::my_rank << ": expected " << count
              << " elements, read " << nread << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    std::exit(EXIT_FAILURE);
  }
  return nread;
}
#endif

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Open(const char* fname, FileMode rw)
//! \brief wrapper for {MPI_File_open} versus {std::fopen} including error check
//! This function must not be called by multiple threads in shared memory parallel regions

int IOWrapper::Open(const char* fname, FileMode rw, bool single_file_per_rank) {
  const char* mode;
  switch (rw) {
    case FileMode::read:
      mode = "rb";
      break;
    case FileMode::write:
      mode = "wb";
      break;
    case FileMode::append:
      mode = "ab";
      break;
    default:
      return false;
  }

#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    int mpi_mode;
    switch (rw) {
      case FileMode::read:
        mpi_mode = MPI_MODE_RDONLY;
        break;
      case FileMode::write:
        mpi_mode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
        MPI_File_delete(fname, MPI_INFO_NULL); // truncation
        break;
      case FileMode::append:
        mpi_mode = MPI_MODE_WRONLY | MPI_MODE_APPEND;
        break;
      default:
        return false;
    }

    int errcode = MPI_File_open(comm_, fname, mpi_mode, MPI_INFO_NULL, &fh_);
    if (errcode != MPI_SUCCESS) {
      char msg[MPI_MAX_ERROR_STRING];
      int resultlen;
      MPI_Error_string(errcode, msg, &resultlen);
      Kokkos::printf("%.*s\n", resultlen, msg);
      MPI_Abort(MPI_COMM_WORLD, 1);
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "File '" << fname << "' could not be opened"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    FILE* local_fh;
    if ((local_fh = std::fopen(fname, mode)) == nullptr) {
      perror("Error opening file");
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "File '" << fname << "' could not be opened"
                << std::endl;
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      std::exit(EXIT_FAILURE);
    }
    fh_ = reinterpret_cast<IOWrapperFile>(local_fh);
  }
#else
  FILE* local_fh;
  if ((local_fh = std::fopen(fname, mode)) == nullptr) {
    perror("Error opening file");
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "File '" << fname << "' could not be opened"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  fh_ = local_fh;
#endif

  return true;
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read_bytes(void *buf, IOWrapperSizeT size, IOWrapperSizeT cnt
//!                                  , bool single_file_per_rank)
//! \brief wrapper for {MPI_File_read} versus {std::fread}.  Returns number of byte-blocks
//! of given "size" actually read.

std::size_t IOWrapper::Read_bytes(void *buf, IOWrapperSizeT size, IOWrapperSizeT cnt,
                                  bool single_file_per_rank, bool require_full_read) {
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    MPI_Status status;
    int errcode = MPI_File_read(fh_, buf, cnt*size, MPI_BYTE, &status);
    if (errcode != MPI_SUCCESS) {
      char msg[MPI_MAX_ERROR_STRING];
      int resultlen;
      MPI_Error_string(errcode, msg, &resultlen);
      Kokkos::printf("%.*s\n", resultlen, msg);
      return 0;
    }
    int nread;
    if (MPI_Get_count(&status,MPI_BYTE,&nread) == MPI_UNDEFINED) {return 0;}
    if (require_full_read && nread != cnt*size) {
      std::cerr << "MPI file read was truncated: expected " << cnt*size
                << " bytes, read " << nread << std::endl;
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      std::exit(EXIT_FAILURE);
    }
    return nread/size;
  } else {
    if (require_full_read) {
      return ReadPerRankFile(buf, size, cnt, reinterpret_cast<FILE*>(fh_));
    }
    return std::fread(buf, size, cnt, reinterpret_cast<FILE*>(fh_));
  }
#else
  return std::fread(buf, size, cnt, fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read_bytes_at(void *buf, IOWrapperSizeT size,
//!                                  IOWrapperSizeT cnt, IOWrapperSizeT offset,
//!                                  bool single_file_per_rank)
//! \brief wrapper for {MPI_File_read_at} versus {std::fseek+std::fread}
//! Returns number of byte-blocks of given "size" actually read.

std::size_t IOWrapper::Read_bytes_at(void *buf, IOWrapperSizeT size,
                                     IOWrapperSizeT cnt, IOWrapperSizeT offset,
                                     bool single_file_per_rank) {
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    MPI_Status status;
    int errcode = MPI_File_read_at(fh_, offset, buf, cnt*size, MPI_BYTE, &status);
    if (errcode != MPI_SUCCESS) {
      char msg[MPI_MAX_ERROR_STRING];
      int resultlen;
      MPI_Error_string(errcode, msg, &resultlen);
      Kokkos::printf("%.*s\n", resultlen, msg);
      return 0;
    }
    int nread;
    if (MPI_Get_count(&status,MPI_BYTE,&nread) == MPI_UNDEFINED) {return 0;}
    return nread/size;
  } else {
    std::fseek(reinterpret_cast<FILE*>(fh_), offset, SEEK_SET);
    return ReadPerRankFile(buf, size, cnt, reinterpret_cast<FILE*>(fh_));
  }
#else
  std::fseek(fh_, offset, SEEK_SET);
  return std::fread(buf, size, cnt, fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read_bytes_at_all(void *buf, IOWrapperSizeT size,
//!                                      IOWrapperSizeT cnt, IOWrapperSizeT offset,
//!                                      bool single_file_per_rank)
//! \brief wrapper for {MPI_File_read_at_all} versus {std::fseek+std::fread}
//! Returns number of byte-blocks of given "size" actually read.

std::size_t IOWrapper::Read_bytes_at_all(void *buf, IOWrapperSizeT size,
                                         IOWrapperSizeT cnt, IOWrapperSizeT offset,
                                         bool single_file_per_rank) {
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    MPI_Status status;
    int errcode = MPI_File_read_at_all(fh_, offset, buf, cnt*size, MPI_BYTE, &status);
    if (errcode != MPI_SUCCESS) {
      char msg[MPI_MAX_ERROR_STRING];
      int resultlen;
      MPI_Error_string(errcode, msg, &resultlen);
      Kokkos::printf("%.*s\n", resultlen, msg);
      return 0;
    }
    int nread;
    if (MPI_Get_count(&status,MPI_BYTE,&nread) == MPI_UNDEFINED) {return 0;}
    return nread/size;
  } else {
    std::fseek(reinterpret_cast<FILE*>(fh_), offset, SEEK_SET);
    return ReadPerRankFile(buf, size, cnt, reinterpret_cast<FILE*>(fh_));
  }
#else
  std::fseek(fh_, offset, SEEK_SET);
  return std::fread(buf, size, cnt, fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read_Reals(void *buf, IOWrapperSizeT cnt,
//!                               bool single_file_per_rank)
//! \brief wrapper for {MPI_File_read} versus {std::fread} for reading Athena Reals.
//! Returns number of Reals actually read.

std::size_t IOWrapper::Read_Reals(void *buf, IOWrapperSizeT cnt,
                                  bool single_file_per_rank) {
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    MPI_Status status;
    int errcode = MPI_File_read(fh_, buf, cnt, MPI_ATHENA_REAL, &status);
    if (errcode != MPI_SUCCESS) {
      char msg[MPI_MAX_ERROR_STRING];
      int resultlen;
      MPI_Error_string(errcode, msg, &resultlen);
      Kokkos::printf("%.*s\n", resultlen, msg);
      return 0;
    }
    int nread;
    if (MPI_Get_count(&status,MPI_ATHENA_REAL,&nread) == MPI_UNDEFINED) {return 0;}
    return nread;
  } else {
    return ReadPerRankFile(buf, sizeof(Real), cnt, reinterpret_cast<FILE*>(fh_));
  }
#else
  return std::fread(buf, sizeof(Real), cnt, fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read_Reals_at(void *buf, IOWrapperSizeT cnt,
//!                                  IOWrapperSizeT offset, bool single_file_per_rank)
//! \brief wrapper for {MPI_File_read_at} versus {std::fseek+std::fread} for reading
//!  Athena Reals in parallel.  Returns number of Reals actually read.

std::size_t IOWrapper::Read_Reals_at(void *buf, IOWrapperSizeT cnt,
                                     IOWrapperSizeT offset, bool single_file_per_rank) {
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    MPI_Status status;
    int errcode = MPI_File_read_at(fh_, offset, buf, cnt, MPI_ATHENA_REAL, &status);
    if (errcode != MPI_SUCCESS) {
      char msg[MPI_MAX_ERROR_STRING];
      int resultlen;
      MPI_Error_string(errcode, msg, &resultlen);
      Kokkos::printf("%.*s\n", resultlen, msg);
      return 0;
    }
    int nread;
    if (MPI_Get_count(&status,MPI_ATHENA_REAL,&nread) == MPI_UNDEFINED) {return 0;}
    return nread;
  } else {
    std::fseek(reinterpret_cast<FILE*>(fh_), offset, SEEK_SET);
    return ReadPerRankFile(buf, sizeof(Real), cnt, reinterpret_cast<FILE*>(fh_));
  }
#else
  std::fseek(fh_, offset, SEEK_SET);
  return std::fread(buf, sizeof(Real), cnt, fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read_Reals_at_all(void *buf, IOWrapperSizeT cnt,
//!                                      IOWrapperSizeT offset, bool single_file_per_rank)
//! \brief wrapper for {MPI_File_read_at_all} versus {std::fseek+std::fread} for reading
//!  Athena Reals in parallel.  Returns number of Reals actually read.

std::size_t IOWrapper::Read_Reals_at_all(void *buf, IOWrapperSizeT cnt,
                                         IOWrapperSizeT offset,
                                         bool single_file_per_rank) {
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    MPI_Status status;
    int errcode = MPI_File_read_at_all(fh_, offset, buf, cnt, MPI_ATHENA_REAL, &status);
    if (errcode != MPI_SUCCESS) {
      char msg[MPI_MAX_ERROR_STRING];
      int resultlen;
      MPI_Error_string(errcode, msg, &resultlen);
      Kokkos::printf("%.*s\n", resultlen, msg);
      return 0;
    }
    int nread;
    if (MPI_Get_count(&status,MPI_ATHENA_REAL,&nread) == MPI_UNDEFINED) {return 0;}
    return nread;
  } else {
    std::fseek(reinterpret_cast<FILE*>(fh_), offset, SEEK_SET);
    return ReadPerRankFile(buf, sizeof(Real), cnt, reinterpret_cast<FILE*>(fh_));
  }
#else
  std::fseek(fh_, offset, SEEK_SET);
  return std::fread(buf, sizeof(Real), cnt, fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Write_any_type()
//! \brief wrapper for {MPI_File_write} versus {std::fwrite} for writing any type of
//! data, specified by the mpitype argument. Returns number of data elements of given
//! "type" actually written.

std::size_t IOWrapper::Write_any_type(const void *buf, IOWrapperSizeT cnt,
                                      std::string datatype, bool single_file_per_rank) {
#if MPI_PARALLEL_ENABLED
  if (single_file_per_rank) {
    return WriteStdio(buf, cnt, datatype, reinterpret_cast<FILE*>(fh_));
  } else {
    const MPI_Datatype mpitype = MPIDataType(datatype);
    MPI_Status status;
    const int error_code = MPI_File_write(fh_, buf, cnt, mpitype, &status);
    if (error_code != MPI_SUCCESS) {
      IOWrapperFatal(MPIErrorMessage("MPI_File_write", error_code));
    }
    return CheckMPIWrite(status, mpitype, cnt, "MPI_File_write");
  }
#else
  return WriteStdio(buf, cnt, datatype, fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Write_any_type_at()
//! \brief wrapper for {MPI_File_write_at} versus {std::fseek+std::fwrite} for writing any
//! type of data, specified by the datatype argument. Returns number of data elements of
//! given "type" actually written.

std::size_t IOWrapper::Write_any_type_at(const void *buf, IOWrapperSizeT cnt,
                                         IOWrapperSizeT offset, std::string datatype,
                                         bool single_file_per_rank) {
#if MPI_PARALLEL_ENABLED
  if (single_file_per_rank) {
    FILE *file = reinterpret_cast<FILE*>(fh_);
    SeekStdioForWrite(file, offset);
    return WriteStdio(buf, cnt, datatype, file);
  } else {
    const MPI_Datatype mpitype = MPIDataType(datatype);
    MPI_Status status;
    const int error_code = MPI_File_write_at(fh_, offset, buf, cnt, mpitype, &status);
    if (error_code != MPI_SUCCESS) {
      IOWrapperFatal(MPIErrorMessage("MPI_File_write_at", error_code));
    }
    return CheckMPIWrite(status, mpitype, cnt, "MPI_File_write_at");
  }
#else
  SeekStdioForWrite(fh_, offset);
  return WriteStdio(buf, cnt, datatype, fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Write_any_type_at_all()
//! \brief wrapper for {MPI_File_write_at_all} versus {std::fseek+std::fwrite} for writing
//! any type of data, specified by the datatype argument.
//! Returns number of data elements of given "type" actually written.

std::size_t IOWrapper::Write_any_type_at_all(const void *buf, IOWrapperSizeT cnt,
                                            IOWrapperSizeT offset, std::string datatype,
                                            bool single_file_per_rank) {
#if MPI_PARALLEL_ENABLED
  if (single_file_per_rank) {
    FILE *file = reinterpret_cast<FILE*>(fh_);
    SeekStdioForWrite(file, offset);
    return WriteStdio(buf, cnt, datatype, file);
  } else {
    const MPI_Datatype mpitype = MPIDataType(datatype);
    MPI_Status status;
    const int error_code = MPI_File_write_at_all(fh_, offset, buf, cnt, mpitype, &status);
    if (error_code != MPI_SUCCESS) {
      IOWrapperFatal(MPIErrorMessage("MPI_File_write_at_all", error_code));
    }
    return CheckMPIWrite(status, mpitype, cnt, "MPI_File_write_at_all");
  }
#else
  SeekStdioForWrite(fh_, offset);
  return WriteStdio(buf, cnt, datatype, fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void IOWrapper::Close(bool single_file_per_rank)
//  \brief wrapper for {MPI_File_close} versus {std::fclose}

int IOWrapper::Close(bool single_file_per_rank) {
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    const int error_code = MPI_File_close(&fh_);
    if (error_code != MPI_SUCCESS) {
      IOWrapperFatal(MPIErrorMessage("MPI_File_close", error_code));
    }
    return 0;
  } else {
    FILE *file = reinterpret_cast<FILE*>(fh_);
    const int close_status = std::fclose(file);
    const int close_errno = errno;
    if (close_status != 0) {
      IOWrapperFatal("could not flush and close output file: " +
                     std::string(std::strerror(close_errno)));
    }
    return 0;
  }
#else
  const int close_status = std::fclose(fh_);
  const int close_errno = errno;
  if (close_status != 0) {
    IOWrapperFatal("could not flush and close output file: " +
                   std::string(std::strerror(close_errno)));
  }
  return 0;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Seek(IOWrapperSizeT offset, bool single_file_per_rank)
//  \brief wrapper for {MPI_File_seek} versus {std::fseek}

int IOWrapper::Seek(IOWrapperSizeT offset, bool single_file_per_rank) {
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    return MPI_File_seek(fh_, offset, MPI_SEEK_SET);
  } else {
    return std::fseek(reinterpret_cast<FILE*>(fh_), offset, SEEK_SET);
  }
#else
  return std::fseek(fh_, offset, SEEK_SET);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn IOWrapperSizeT IOWrapper::GetPosition(bool single_file_per_rank)
//  \brief wrapper for {MPI_File_get_position} versus {ftell}

IOWrapperSizeT IOWrapper::GetPosition(bool single_file_per_rank) {
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    MPI_Offset position;
    MPI_File_get_position(fh_, &position);
    return position;
  } else {
    int64_t pos = ftell(reinterpret_cast<FILE*>(fh_));
    return pos;
  }
#else
  int64_t pos = ftell(fh_);
  return pos;
#endif
}
