#!/bin/bash

# AMD Ryzen 9950X (Zen5) + GCC 15.2.0 + OpenMPI 安全构建脚本
# 基于 AthenaK 构建文档: https://github.com/IAS-Astrophysics/athenak/wiki/Build

# 创建构建目录
mkdir -p build_9950x
cd build_9950x

# 检查编译器版本和支持的指令集
echo "检查编译器版本..."
mpicc --version | head -1
mpicxx --version | head -1

# 测试编译器是否支持特定指令集
echo "测试编译器支持的指令集..."

# 测试 -mprefetchwt1 是否支持
if mpicc -mprefetchwt1 -x c - -o /dev/null 2>/dev/null <<< 'int main(){return 0;}'; then
    PREFETCH_FLAG="-mprefetchwt1"
    echo "✓ 支持 -mprefetchwt1"
else
    PREFETCH_FLAG=""
    echo "✗ 不支持 -mprefetchwt1，将跳过此选项"
fi

# 测试 -mavx512vnni 是否支持
if mpicc -mavx512vnni -x c - -o /dev/null 2>/dev/null <<< 'int main(){return 0;}'; then
    AVX512VNNI_FLAG="-mavx512vnni"
    echo "✓ 支持 -mavx512vnni"
else
    AVX512VNNI_FLAG=""
    echo "✗ 不支持 -mavx512vnni，将跳过此选项"
fi

# 构建编译器标志
CXX_FLAGS="-march=znver5 -O3 -ffast-math -funroll-loops -flto"
C_FLAGS="-march=znver5 -O3 -ffast-math -funroll-loops -flto"

if [ -n "$PREFETCH_FLAG" ]; then
    CXX_FLAGS="$CXX_FLAGS $PREFETCH_FLAG"
    C_FLAGS="$C_FLAGS $PREFETCH_FLAG"
fi

if [ -n "$AVX512VNNI_FLAG" ]; then
    CXX_FLAGS="$CXX_FLAGS $AVX512VNNI_FLAG"
    C_FLAGS="$C_FLAGS $AVX512VNNI_FLAG"
fi

echo "最终使用的编译器标志:"
echo "CXX_FLAGS: $CXX_FLAGS"
echo "C_FLAGS: $C_FLAGS"

# 针对 AMD Ryzen 9950X (Zen5) 优化的 cmake 配置
cmake \
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D CMAKE_C_COMPILER=mpicc \
    -D Kokkos_ARCH_NATIVE=ON \
    -D Athena_ENABLE_MPI=ON \
    -D CMAKE_BUILD_TYPE=Release \
    -D CMAKE_CXX_FLAGS="$CXX_FLAGS" \
    -D CMAKE_C_FLAGS="$C_FLAGS" \
    -D CMAKE_EXE_LINKER_FLAGS="-flto" \
    -D Kokkos_ENABLE_OPENMP=ON \
    -D Kokkos_ENABLE_SERIAL=ON \
    -D Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION=ON \
    -D Kokkos_ENABLE_DEPRECATED_CODE_4=OFF \
    -D Kokkos_ENABLE_DEBUG=OFF \
    -D Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=OFF \
    -D Kokkos_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK=OFF \
    -D Kokkos_ENABLE_TUNING=OFF \
    -D PROBLEM=gr_torus \
    -D Athena_SINGLE_PRECISION=ON \
    ../

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ cmake 配置成功！"
    echo "针对 AMD Ryzen 9950X (Zen5) 架构优化"
    echo "启用的指令集: 通过 -march=znver5 自动检测 Zen5 架构支持的所有指令集"
    echo "现在可以运行: make -j$(nproc)"
else
    echo ""
    echo "✗ cmake 配置失败！"
    echo "请检查错误信息并修复问题"
    exit 1
fi
