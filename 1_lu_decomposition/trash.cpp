/*

template<class T>
using Blocks = std::vector<std::vector<T>>;

template<class T>
Blocks<T**> SplitBlocks(T** matrix, size_t n, size_t blocks) {

    const size_t blockSize = n / blocks;

    Blocks<T**> ret;
    ret.resize(blocks);
    for (auto& block : ret) {
        block.resize(blocks);
        for (auto &matrix : block) {
            matrix = Allocate<T>(blockSize);
        }
    }

    if (matrix) {
        for (size_t i = 0; i < blocks; ++i) {
            for (size_t j = 0; j < blocks; ++j) {
                Copy<T>(ret[i][j], matrix, blockSize, j * blockSize, i * blockSize);
            }
        }
    } else {
        for (size_t i = 0; i < blocks; ++i) {
            for (size_t j = 0; j < blocks; ++j) {
                FillZeros<T>(ret[i][j], blockSize);
            }
        }
    }

    return ret;
}

template<class T>
void GlueBlocks(T** ret, Blocks<T**> blocks, size_t n, size_t blocksSize) {
    const auto blockSize = n / blocksSize;
    for (size_t i = 0; i < blocksSize; ++i) {
        for (size_t j = 0; j < blocksSize; ++j) {
            Copy(ret, blocks[i][j], blockSize, 0, 0, j * blockSize, i * blockSize);
        }
    }
}

template<class T>
void Decompose(T** matrix, size_t size) {
    for (size_t i = 1; i < size; ++i) {
        for (size_t k = 0; k <= i - 1; ++k) {
            matrix[i][k] = matrix[i][k] / matrix[k][k];
            for (size_t j = k + 1; j < size; ++j) {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
        }
    }
}

template<class T>
void MP_Decompose_Block(T** matrix, T** L, T** U, size_t size) {
    const auto blocksSize = CalculateBlocksSize(size);
    // std::cout << "Blocks size: " << blocksSize << std::endl;
    if(blocksSize == 1) {
        Decompose(matrix, size);
        return;
    }

    auto ABlocks = SplitBlocks(matrix, size, blocksSize);
    auto LBlocks = SplitBlocks<double>(nullptr, size, blocksSize);
    auto UBlocks = SplitBlocks<double>(nullptr, size, blocksSize);

    std::vector<std::thread> procs;

    _ASSERT(blocksSize == 2 || blocksSize == 1);

    const size_t blockSize = size / blocksSize;
    // std::cout << "Block size: " << blockSize << std::endl;
    for (size_t i = 0; i < blocksSize; ++i) {
        for (size_t j = 0; j < blocksSize; ++j) {
            procs.emplace_back([&, i, j]() {
                // std::cout << "A" << j << i << std::endl;
                // Print(ABlocks[i][j], blockSize);

                Decompose(ABlocks[i][j], blockSize);
                auto res = ExtractDecomposition(ABlocks[i][j], blockSize);
                if(j != 1 || i != 0)
                {
                    LBlocks[i][j] = std::get<0>(res);
                }
                if(j != 0 || i != 1)
                {
                    UBlocks[i][j] = std::get<1>(res);
                }

                // std::cout << "L" << j << i << std::endl;
                // Print(LBlocks[i][j], blockSize);
                // std::cout << "U" << j << i << std::endl;
                // Print(UBlocks[i][j], blockSize);

            });
            // procs.back().join();
        }
    }

    for (auto& proc: procs) {
        proc.join();
    }

    GlueBlocks(L, LBlocks, size, blocksSize);
    GlueBlocks(U, UBlocks, size, blocksSize);
}

void LU_Decomposition(double* A, double* L, double* U, int n)
{
    auto AConv = Convert(A, n);
    auto LConv = Convert(L, n);
    auto UConv = Convert(U, n);

    MP_Decompose_Block(AConv, LConv, UConv, n);

    DeallocInterpret(AConv);
    Copy(L, LConv, n);
    Copy(U, UConv, n);

    Dealloc(LConv, n);
    Dealloc(UConv, n);
}

*/