using System;

namespace Math.Core
{
    public partial class Functions
    {
        /// <summary>
        /// Cholesky 分解：A = L * L^T，返回 0 表示成功，1 表示矩阵非正定。
        /// </summary>
        /// <param name="input">对称正定矩阵 A（只使用上三角）。</param>
        /// <param name="output">下三角矩阵 L（仅下三角会被写入）。</param>
        public static int Cholesky(double[][] input, double[][] output)
        {
            if (input == null || output == null)
            {
                throw new Exception("matrix is null");
            }

            int size = input.Length;
            if (size == 0 || output.Length != size)
            {
                throw new Exception("matrix size mismatch");
            }

            for (int i = 0; i < size; i++)
            {
                if (input[i].Length != size || output[i].Length != size)
                {
                    throw new Exception("matrix must be square");
                }
            }

            for (int i = 0; i < size; i++)
            {
                for (int j = i; j < size; j++)
                {
                    double sum = input[i][j];
                    for (int k = i - 1; k >= 0; k--)
                    {
                        sum -= output[i][k] * output[j][k];
                    }

                    if (i == j)
                    {
                        if (sum <= 0.0)
                        {
                            return 1;
                        }

                        output[i][i] = Math.Sqrt(sum);
                    }
                    else
                    {
                        output[j][i] = sum / output[i][i];
                    }
                }
            }

            return 0;
        }
    }
}
