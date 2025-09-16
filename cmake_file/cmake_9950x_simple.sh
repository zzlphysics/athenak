#!/bin/bash

# AMD Ryzen 9950X (Zen5) + GCC 15.2.0 + OpenMPI 简化构建脚本
# 基于 AthenaK 构建文档: https://github.com/IAS-Astrophysics/athenak/wiki/Build

# 创建构建目录
mkdir -p build_9950x
cd build_9950x

# 针对 AMD Ryzen 9950X (Zen5) 优化的 cmake 配置
cmake \
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D CMAKE_C_COMPILER=mpicc \
    -D Kokkos_ARCH_ZEN5=ON \
    -D Kokkos_ARCH_NATIVE=ON \
    -D Athena_ENABLE_MPI=ON \
    -D CMAKE_BUILD_TYPE=Release \
    -D CMAKE_CXX_FLAGS="-march=znver5 -O3 -ffast-math -funroll-loops -flto" \
    -D CMAKE_C_FLAGS="-march=znver5 -O3 -ffast-math -funroll-loops -flto" \
    -D CMAKE_EXE_LINKER_FLAGS="-flto" \
    -D Kokkos_ENABLE_OPENMP=ON \
    -D Kokkos_ENABLE_SERIAL=ON \
    -D Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION=ON \
    -D PROBLEM=gr_torus \
    ../

echo "cmake 配置完成！"
echo "针对 AMD Ryzen 9950X (Zen5) 架构优化"
echo "启用的指令集: 通过 -march=znver5 自动检测 Zen5 架构支持的所有指令集"
echo "现在可以运行: make -j$(nproc)"