using System;
using Math.Core.Extensions;
using MathNet.Numerics.LinearAlgebra;

namespace Math.Core
{
    public partial class Functions
    {
        /// <summary>
        /// QR分解，将输入矩阵A分解为一个正交矩阵Q和一个上三角矩阵R
        /// </summary>
        /// <param name="A"></param>
        /// <param name="Q"></param>
        /// <param name="R"></param>
        /// <param name="perm"></param>
        /// <param name="sort"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        public static int QRfactorize(double[][] A, out double[][] Q, out double[][] R, out int[] perm, int sort = 0, double eps = 1.00e-16)
        {
            int rows = A.Length;
            int cols = A[0].Length;
            perm = new int[cols];
            for (int i = 0; i < cols; i++)
            {
                perm[i] = i;
            }

            if (sort == 0 && rows > 0 && cols > 0)
            {
                var matrix = A.ToMatrix();
                var qr = matrix.QR();
                Q = ToJagged(qr.Q);
                R = ToJagged(qr.R);
                return ComputeRank(qr.R, eps);
            }

            if (rows == 0)
            {
                Q = new double[0][].Resize(0, 0);
                R = new double[0][].Resize(0, cols);
                return 0;
            }

            if (rows == 1)
            {
                R = new double[0][].Resize(rows, cols);
                for (int j = 0; j < cols; j++)
                {
                    R[0][j] = A[0][j];
                }

                double sum = 0.0;
                for (int j = 0; j < cols; j++)
                {
                    double d = A[0][j];
                    sum += d * d;
                }

                Q = new double[0][].Resize(rows, rows);
                Q[0][0] = 1.0;
                return sum > 0.0 ? 1 : 0;
            }

            R = new double[0][].Resize(rows, cols);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    R[i][j] = 0.0;
                }
            }

            var C = new double[0][].Resize(rows, cols);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    C[i][j] = A[i][j];
                }
            }

            int rank = QRfactorize(ref C, out var diag, out perm, sort, eps);
            QRcalcHouseholders(out Q, C, perm, rank);
            QRcalcUpperRightTriangle(out R, C, diag, perm);

            return rank;
        }

        /// <summary>
        /// QR分解算法，用于将矩阵C分解为一个正交矩阵Q和一个上三角矩阵R
        /// </summary>
        /// <param name="C"></param>
        /// <param name="diagR"></param>
        /// <param name="perm"></param>
        /// <param name="sort"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        public static int QRfactorize(ref double[][] C, out double[] diagR, out int[] perm, int sort = 0, double eps = 1.00e-16)
        {
            int rows = C.Length;
            int cols = C[0].Length;
            int minDim = Math.Min(rows, cols);

            int rank = 0;
            diagR = new double[minDim];

            double refElC = 1.0;
            if (rows > 0)
            {
                refElC = (Math.Abs(C[0][0]) + Math.Abs(C[0][cols - 1]) + Math.Abs(C[rows / 2][cols / 2])
                    + Math.Abs(C[rows - 1][0]) + Math.Abs(C[rows - 1][cols - 1])) / 5.0;
            }

            perm = new int[cols];
            for (int i = 0; i < cols; i++)
            {
                perm[i] = i;
            }

            for (int l = 0; l < minDim; l++)
            {
                int pivot = perm[l];
                int jMax = l;
                double gamma = 0.0;
                for (int i = l; i < rows; i++)
                {
                    double c = C[i][pivot];
                    gamma += c * c;
                }

                gamma = Math.Sqrt(gamma);

                if (sort != 0)
                {
                    for (int j = l + 1; j < cols; j++)
                    {
                        int pj = perm[j];
                        double sum = 0.0;
                        for (int i = l; i < rows; i++)
                        {
                            double c = C[i][pj];
                            sum += c * c;
                        }

                        sum = Math.Sqrt(sum);
                        if (sum > gamma)
                        {
                            gamma = sum;
                            jMax = j;
                        }
                    }

                    perm[l] = perm[jMax];
                    perm[jMax] = pivot;
                    pivot = perm[l];
                }

                if (gamma > eps * refElC * Math.Max(rows, cols))
                {
                    rank++;
                    double d = C[l][pivot];
                    gamma *= (d < 0.0) ? -1.0 : 1.0;
                    double beta = Math.Sqrt(2 * gamma * (gamma + d));
                    diagR[l] = -gamma;

                    C[l][pivot] += gamma;
                    for (int k = l; k < rows; k++)
                    {
                        C[k][pivot] /= beta;
                    }

                    for (int j = l + 1; j < cols; j++)
                    {
                        int pj = perm[j];
                        double p = 0.0;
                        for (int k = l; k < rows; k++)
                        {
                            p += C[k][pivot] * C[k][pj];
                        }

                        for (int i = l; i < rows; i++)
                        {
                            C[i][pj] -= 2 * p * C[i][pivot];
                        }
                    }
                }
                else
                {
                    diagR[l] = 0.0;
                    for (int k = l; k < rows; k++)
                    {
                        C[k][pivot] = 0.0;
                    }

                    if (sort != 0)
                    {
                        var trimmed = new double[rank];
                        Array.Copy(diagR, trimmed, rank);
                        diagR = trimmed;
                        break;
                    }
                }
            }

            return rank;
        }

        /// <summary>
        /// QR分解中的Householder变换
        /// </summary>
        /// <param name="Q"></param>
        /// <param name="C"></param>
        /// <param name="perm"></param>
        /// <param name="rank"></param>
        /// <exception cref="Exception"></exception>
        public static void QRcalcHouseholders(out double[][] Q, double[][] C, int[] perm, int rank)
        {
            int rows = C.Length;
            int cols = C[0].Length;
            int minDim = Math.Min(rows, cols);
            if (minDim < rank)
            {
                throw new Exception("wrong rank");
            }

            Q = new double[0][].Resize(rows, rows);
            unit_matrix(ref Q);

            for (int l = 0; l < rank; l++)
            {
                int pivot = perm[l];
                for (int j = 0; j < rows; j++)
                {
                    double p = 0.0;
                    for (int k = l; k < rows; k++)
                    {
                        p += C[k][pivot] * Q[j][k];
                    }

                    for (int k = l; k < rows; k++)
                    {
                        Q[j][k] -= 2 * p * C[k][pivot];
                    }
                }
            }
        }

        /// <summary>
        /// 通过QR分解法计算上三角矩阵
        /// </summary>
        /// <param name="R"></param>
        /// <param name="C"></param>
        /// <param name="diagR"></param>
        /// <param name="perm"></param>
        public static void QRcalcUpperRightTriangle(out double[][] R, double[][] C, double[] diagR, int[] perm)
        {
            int rows = C.Length;
            int cols = C[0].Length;
            int minDim = Math.Min(rows, cols);
            int ndia = diagR.Length;

            R = new double[0][].Resize(rows, cols);
            for (int i = 0; i < minDim; i++)
            {
                R[i][i] = (i < ndia) ? diagR[i] : 0.0;
                for (int j = i + 1; j < cols; j++)
                {
                    R[i][j] = C[i][perm[j]];
                }
            }
        }

        /// <summary>
        /// 利用QR分解的结果，通过Householder变换将输入向量 x 转换为新的向量 y
        /// </summary>
        /// <param name="y"></param>
        /// <param name="C"></param>
        /// <param name="perm"></param>
        /// <param name="rank"></param>
        /// <param name="x"></param>
        /// <exception cref="Exception"></exception>
        public static void QRapplyHouseholders(out double[] y, double[][] C, int[] perm, int rank, double[] x)
        {
            int rows = C.Length;
            if (x.Length != rows)
            {
                throw new Exception("x and Q are incompatible");
            }

            int minDim = Math.Min(rows, C[0].Length);
            if (minDim < rank)
            {
                throw new Exception("wrong rank");
            }

            y = (double[])x.Clone();
            for (int l = 0; l < rank; l++)
            {
                int pivot = perm[l];
                double p = 0.0;
                for (int k = l; k < rows; k++)
                {
                    p += C[k][pivot] * y[k];
                }

                for (int i = l; i < rows; i++)
                {
                    y[i] -= 2 * p * C[i][pivot];
                }
            }
        }

        /// <summary>
        /// 利用QR分解的结果，通过回代法计算线性方程组的解向量 x
        /// </summary>
        /// <param name="x"></param>
        /// <param name="C"></param>
        /// <param name="diagR"></param>
        /// <param name="perm"></param>
        /// <param name="rank"></param>
        /// <param name="rhs"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static int QRbacksubstitute(out double[] x, double[][] C, double[] diagR, int[] perm, int rank, double[] rhs)
        {
            int rows = C.Length;
            int cols = C[0].Length;
            int minDim = Math.Min(rows, cols);

            if (minDim < rank)
            {
                throw new Exception("wrong rank");
            }

            if (rhs.Length != rows)
            {
                throw new Exception("right hand side and matrix are incompatible");
            }

            x = new double[cols];
            for (int i = rank - 1; i >= 0; i--)
            {
                double sum = 0.0;
                for (int k = i + 1; k < rank; k++)
                {
                    sum += C[i][perm[k]] * x[perm[k]];
                }

                double rii = diagR[i];
                x[perm[i]] = (rhs[i] - sum) / rii;
            }

            return minDim - rank;
        }

        /// <summary>
        /// 利用QR分解的结果，通过回代法计算线性方程组的解向量 x
        /// </summary>
        /// <param name="x"></param>
        /// <param name="R"></param>
        /// <param name="perm"></param>
        /// <param name="rank"></param>
        /// <param name="rhs"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static int QRbacksubstitute(out double[] x, double[][] R, int[] perm, int rank, double[] rhs)
        {
            int rows = R.Length;
            int cols = R[0].Length;
            int minDim = Math.Min(rows, cols);

            if (minDim < rank)
            {
                throw new Exception("wrong rank");
            }

            if (rhs.Length != rows)
            {
                throw new Exception("right hand side and matrix are incompatible");
            }

            x = new double[cols];
            for (int i = rank - 1; i >= 0; i--)
            {
                double sum = 0.0;
                for (int k = i + 1; k < rank; k++)
                {
                    sum += R[i][k] * x[perm[k]];
                }

                double rii = R[i][i];
                x[perm[i]] = (rhs[i] - sum) / rii;
            }

            return minDim - rank;
        }

        /// <summary>
        /// 使用 QR 分解将矩阵 A 分解为一个正交矩阵 Q 和一个上三角矩阵 R，然后通过回代法求解线性方程组
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <param name="Q"></param>
        /// <param name="R"></param>
        /// <param name="perm"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static int QRsolve(double[][] A, double[] b, ref double[] x, ref double[][] Q, ref double[][] R, out int[] perm, double eps = 1.00e-16)
        {
            int rows = A.Length;
            int cols = A[0].Length;
            int minDim = Math.Min(rows, cols);

            if (b.Length != rows)
            {
                throw new Exception("size(b) != rows(A)");
            }

            var matrix = A.ToMatrix();
            var qr = matrix.QR();
            Q = ToJagged(qr.Q);
            R = ToJagged(qr.R);
            perm = new int[cols];
            for (int i = 0; i < cols; i++)
            {
                perm[i] = i;
            }

            x = qr.Solve(Vector<double>.Build.DenseOfArray(b)).ToArray();
            int rank = ComputeRank(qr.R, eps);
            return minDim - rank;
        }

        /// <summary>
        /// 通过 QR 分解求解线性方程组
        /// </summary>
        /// <param name="x"></param>
        /// <param name="perm"></param>
        /// <param name="diagR"></param>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static int QRsolveWithDiagR(out double[] x, out int[] perm, out double[] diagR, ref double[][] A, double[] b, double eps = 1.00e-16)
        {
            int rows = A.Length;
            int cols = A[0].Length;
            int minDim = Math.Min(rows, cols);

            if (b.Length != rows)
            {
                throw new Exception("size(b) != rows(A)");
            }

            x = new double[cols];
            int rank = QRfactorize(ref A, out diagR, out perm, 1, eps);
            if (rank == 0)
            {
                return minDim;
            }

            QRapplyHouseholders(out var rhs, A, perm, rank, b);
            QRbacksubstitute(out x, A, diagR, perm, rank, rhs);

            return minDim - rank;
        }

        /// <summary>
        /// 利用QR分解方法来求解线性方程组
        /// </summary>
        /// <param name="x"></param>
        /// <param name="perm"></param>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        public static int QRsolve(ref double[] x, out int[] perm, ref double[][] A, double[] b, double eps = 1.00e-16)
        {
            return QRsolveWithDiagR(out x, out perm, out _, ref A, b, eps);
        }

        /// <summary>
        /// 计算Givens旋转矩阵的参数，即cos(c) 和 sin(s) 值。这些参数用于将一个二维向量旋转到某个特定方向
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <param name="s"></param>
        /// <param name="eps"></param>
        public static void getRotation(double a, double b, out double c, out double s, double eps = 1.00e-12)
        {
            if ((Math.Abs(b) < eps) && (Math.Abs(a) < eps))
            {
                c = 1.0;
                s = 0.0;
            }
            else if (Math.Abs(b) > Math.Abs(a))
            {
                double tau = a / b;
                if (Math.Abs(tau) < eps)
                {
                    s = 1.0;
                    c = 0.0;
                }
                else
                {
                    s = 1.0 / Math.Sqrt(1.0 + tau * tau);
                    c = s * tau;
                }
            }
            else
            {
                double tau = b / a;
                if (Math.Abs(tau) < eps)
                {
                    c = 1.0;
                    s = 0.0;
                }
                else
                {
                    c = 1.0 / Math.Sqrt(1.0 + tau * tau);
                    s = c * tau;
                }
            }
        }

        /// <summary>
        /// 在QR分解过程中删除矩阵R的第k列，并相应地更新正交矩阵Q
        /// </summary>
        /// <param name="Q"></param>
        /// <param name="R"></param>
        /// <param name="k"></param>
        /// <param name="eps"></param>
        /// <exception cref="Exception"></exception>
        public static void QRdelete(ref double[][] Q, ref double[][] R, int k, double eps = 1.00e-16)
        {
            int rows = R.Length;
            int cols = R[0].Length;

            int qRows = Q.Length;
            int qCols = Q[0].Length;
            if (qRows != qCols)
            {
                throw new Exception("cols(Q) != rows(Q)");
            }

            if (rows != qRows)
            {
                throw new Exception("dimensions of Q and R don't match");
            }

            if (k >= cols)
            {
                throw new Exception("k >= cols(R)");
            }

            for (int j = k; j < cols - 1; j++)
            {
                for (int i = 0; i < rows; i++)
                {
                    R[i][j] = R[i][j + 1];
                }
            }

            R = new double[0][].Resize(rows, cols - 1);

            int min = Math.Min(rows - 1, cols - 1);
            for (int j = k; j < min; j++)
            {
                double a = R[j][j];
                double b = R[j + 1][j];
                double r2 = a * a + b * b;
                if (r2 <= 2.0 * eps * eps)
                {
                    continue;
                }

                getRotation(a, b, out double c, out double s, eps);

                for (int i = 0; i < cols - 1; i++)
                {
                    double d = c * R[j][i] + s * R[j + 1][i];
                    double e = c * R[j + 1][i] - s * R[j][i];
                    R[j][i] = d;
                    R[j + 1][i] = e;
                }

                for (int i = 0; i < rows; i++)
                {
                    double d = c * Q[i][j] + s * Q[i][j + 1];
                    double e = c * Q[i][j + 1] - s * Q[i][j];
                    Q[i][j] = d;
                    Q[i][j + 1] = e;
                }
            }
        }

        /// <summary>
        /// 在QR分解过程中将一个新的列向量添加到矩阵R的右侧，并相应地更新正交矩阵Q
        /// </summary>
        /// <param name="Q"></param>
        /// <param name="R"></param>
        /// <param name="col"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static int QRappend(ref double[][] Q, ref double[][] R, double[] col, double eps = 1.00e-16)
        {
            int rows = R.Length;
            int cols = R[0].Length;

            int qRows = Q.Length;
            int qCols = Q[0].Length;
            int colRows = col.Length;
            if (qRows != qCols)
            {
                throw new Exception("cols(Q) != rows(Q)");
            }

            if (rows != qRows)
            {
                throw new Exception("dimensions of Q and R don't match");
            }

            if (colRows != qRows)
            {
                throw new Exception("dimensions of Q and col don't match");
            }

            int ret = 0;
            R = new double[0][].Resize(rows, cols + 1);

            for (int i = 0; i < rows; i++)
            {
                double dummy = 0;
                for (int j = 0; j < rows; j++)
                {
                    dummy += Q[j][i] * col[j];
                }

                R[i][cols] = dummy;
            }

            if (cols >= rows)
            {
                return ret;
            }

            double gamma = 0.0;
            for (int i = cols; i < rows; i++)
            {
                double el = R[i][cols];
                gamma += el * el;
            }

            gamma = Math.Sqrt(gamma);
            if (gamma > eps * Math.Abs(R[0][0]) * Math.Max(cols, rows))
            {
                ret += 1;

                double diag = R[cols][cols];
                gamma *= (diag < 0.0) ? -1.0 : 1.0;

                double beta = Math.Sqrt(2.0 * gamma * (gamma + diag));
                R[cols][cols] += gamma;
                for (int i = cols; i < rows; i++)
                {
                    R[i][cols] /= beta;
                }

                for (int j = 0; j < rows; j++)
                {
                    double p = 0.0;
                    for (int i = cols; i < rows; i++)
                    {
                        p += R[i][cols] * Q[j][i];
                    }

                    for (int i = cols; i < rows; i++)
                    {
                        Q[j][i] -= 2 * p * R[i][cols];
                    }
                }

                R[cols][cols] = -gamma;
            }
            else
            {
                R[cols][cols] = 0.0;
            }

            for (int i = cols + 1; i < rows; i++)
            {
                R[i][cols] = 0.0;
            }

            return ret;
        }

        /// <summary>
        /// 在QR分解过程中将一个新的行向量添加到矩阵R的顶部，并相应地更新正交矩阵Q
        /// </summary>
        /// <param name="Q"></param>
        /// <param name="R"></param>
        /// <param name="row"></param>
        /// <param name="eps"></param>
        /// <exception cref="Exception"></exception>
        public static void QRprependRow(ref double[][] Q, double[][] R, double[] row, double eps = 1.00e-16)
        {
            int rows = R.Length;
            int cols = R[0].Length;

            int qRows = Q.Length;
            int qCols = Q[0].Length;
            int rowCols = row.Length;

            bool updateQ = qRows * qCols > 0;
            if (updateQ)
            {
                if (qRows != qCols)
                {
                    throw new Exception("cols(Q) != rows(Q)");
                }

                if (rows - 1 != qRows)
                {
                    throw new Exception("dimensions of Q and R don't match");
                }
            }

            if (cols != rowCols)
            {
                throw new Exception("dimensions of R and row don't match");
            }

            int originalRows = rows - 1;
            if (originalRows < 0)
            {
                throw new Exception("R must be preallocated with one extra row");
            }

            var rowCopy = (double[])row.Clone();

            if (updateQ)
            {
                Q = new double[0][].Resize(rows, rows);
                for (int i = originalRows; i > 0; i--)
                {
                    for (int j = originalRows; j > 0; j--)
                    {
                        Q[i][j] = Q[i - 1][j - 1];
                    }

                    Q[i][0] = 0.0;
                }

                for (int j = originalRows; j > 0; j--)
                {
                    Q[0][j] = 0.0;
                }

                Q[0][0] = 1.0;
            }

            for (int i = originalRows - 1; i >= 0; i--)
            {
                for (int j = 0; j < cols; j++)
                {
                    R[i + 1][j] = R[i][j];
                }
            }

            int min = Math.Min(originalRows, cols);
            for (int j = 0; j < min; j++)
            {
                double a = rowCopy[j];
                double b = R[j][j];
                double r2 = a * a + b * b;
                if (r2 <= eps * eps)
                {
                    continue;
                }

                getRotation(a, b, out double c, out double s, eps);
                for (int i = 0; i < cols; i++)
                {
                    double d = c * rowCopy[i] + s * R[j][i];
                    double e = c * R[j][i] - s * rowCopy[i];
                    R[j][i] = d;
                    rowCopy[i] = e;
                }

                if (updateQ)
                {
                    for (int i = 0; i < rows; i++)
                    {
                        double d = c * Q[i][j] + s * Q[i][j + 1];
                        double e = c * Q[i][j + 1] - s * Q[i][j];
                        Q[i][j] = d;
                        Q[i][j + 1] = e;
                    }
                }
            }
        }

        /// <summary>
        /// QR分解过程中替换矩阵R的某一列，并相应地更新正交矩阵Q
        /// </summary>
        /// <param name="Q"></param>
        /// <param name="R"></param>
        /// <param name="col"></param>
        /// <param name="k"></param>
        /// <param name="eps"></param>
        /// <exception cref="Exception"></exception>
        public static void QRreplace(ref double[][] Q, ref double[][] R, double[] col, int k, double eps = 1.00e-16)
        {
            int rows = R.Length;
            int cols = R[0].Length;

            int qRows = Q.Length;
            int qCols = Q[0].Length;
            int colRows = col.Length;

            if (qRows != qCols)
            {
                throw new Exception("cols(Q) != rows(Q)");
            }

            if (rows != qRows)
            {
                throw new Exception("dimensions of Q and R don't match");
            }

            if (rows != colRows)
            {
                throw new Exception("dimensions of R and col don't match");
            }

            if (k >= cols)
            {
                throw new Exception("column to be replaced does not exist");
            }

            for (int i = 0; i < rows; i++)
            {
                double dummy = 0;
                for (int j = 0; j < rows; j++)
                {
                    dummy += Q[j][i] * col[j];
                }

                R[i][k] = dummy;
            }

            for (int i = rows - 1; i > k; i--)
            {
                double a = R[i - 1][k];
                double b = R[i][k];
                double r2 = a * a + b * b;
                if (r2 <= eps * eps)
                {
                    continue;
                }

                getRotation(a, b, out double c, out double s, eps);
                for (int j = 0; j < cols; j++)
                {
                    double d = c * R[i - 1][j] + s * R[i][j];
                    double e = c * R[i][j] - s * R[i - 1][j];
                    R[i - 1][j] = d;
                    R[i][j] = e;
                }

                for (int j = 0; j < rows; j++)
                {
                    double d = c * Q[j][i - 1] + s * Q[j][i];
                    double e = c * Q[j][i] - s * Q[j][i - 1];
                    Q[j][i - 1] = d;
                    Q[j][i] = e;
                }
            }

            int min = Math.Min(rows - 1, cols);
            for (int j = k + 1; j < min; j++)
            {
                double a = R[j][j];
                double b = R[j + 1][j];
                double r2 = a * a + b * b;
                if (r2 <= eps * eps)
                {
                    continue;
                }

                getRotation(a, b, out double c, out double s, eps);
                for (int i = 0; i < cols; i++)
                {
                    double d = c * R[j][i] + s * R[j + 1][i];
                    double e = c * R[j + 1][i] - s * R[j][i];
                    R[j][i] = d;
                    R[j + 1][i] = e;
                }

                for (int i = 0; i < rows; i++)
                {
                    double d = c * Q[i][j] + s * Q[i][j + 1];
                    double e = c * Q[i][j + 1] - s * Q[i][j];
                    Q[i][j] = d;
                    Q[i][j + 1] = e;
                }
            }
        }

        private static double[][] ToJagged(Matrix<double> matrix)
        {
            var result = new double[matrix.RowCount][];
            for (int i = 0; i < matrix.RowCount; i++)
            {
                result[i] = matrix.Row(i).ToArray();
            }

            return result;
        }

        private static int ComputeRank(Matrix<double> r, double eps)
        {
            int n = Math.Min(r.RowCount, r.ColumnCount);
            double refEl = n > 0 ? Math.Abs(r[0, 0]) : 1.0;
            double threshold = eps * refEl * Math.Max(r.RowCount, r.ColumnCount);
            int rank = 0;
            for (int i = 0; i < n; i++)
            {
                if (Math.Abs(r[i, i]) > threshold)
                {
                    rank++;
                }
            }

            return rank;
        }
    }
}
