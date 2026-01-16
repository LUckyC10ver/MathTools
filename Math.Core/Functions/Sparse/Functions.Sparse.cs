using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Math.Core
{
    public partial class Functions
    {
        /// <summary>
        /// BiCGSTAB 求解稀疏线性方程组。
        /// </summary>
        public static int BiCGSTAB(
            SparseMatrix A,
            double[] x,
            double[] b,
            Func<Vector<double>, Vector<double>> preconditioner,
            ref int maxIter,
            ref double tol)
        {
            if (A == null)
            {
                throw new Exception("matrix is null");
            }

            if (x == null || b == null)
            {
                throw new Exception("vector is null");
            }

            int n = x.Length;
            if (b.Length != n || A.RowCount != n || A.ColumnCount != n)
            {
                throw new Exception("dimensions of A, x, b do not match");
            }

            var vectorBuilder = Vector<double>.Build;
            var xVec = vectorBuilder.DenseOfArray(x);
            var bVec = vectorBuilder.DenseOfArray(b);
            var r = bVec - A * xVec;
            var rtilde = r.Clone();

            double normb = bVec.L2Norm();
            if (normb == 0.0)
            {
                normb = 1.0;
            }

            double resid = r.L2Norm() / normb;
            if (resid <= tol)
            {
                tol = resid;
                maxIter = 0;
                CopyVector(xVec, x);
                return 0;
            }

            double rho1 = 0.0;
            double rho2 = 1.0;
            double alpha = 1.0;
            double omega = 1.0;

            var p = vectorBuilder.Dense(n, 0.0);
            var v = vectorBuilder.Dense(n, 0.0);

            for (int i = 1; i <= maxIter; i++)
            {
                rho1 = rtilde.DotProduct(r);
                if (rho1 == 0.0)
                {
                    tol = r.L2Norm() / normb;
                    return 2;
                }

                if (i == 1)
                {
                    p = r.Clone();
                }
                else
                {
                    double beta = (rho1 / rho2) * (alpha / omega);
                    p = r + beta * (p - omega * v);
                }

                var phat = ApplyPreconditioner(preconditioner, p);
                v = A * phat;
                alpha = rho1 / rtilde.DotProduct(v);
                var s = r - alpha * v;

                resid = s.L2Norm() / normb;
                if (resid < tol)
                {
                    xVec += alpha * phat;
                    tol = resid;
                    maxIter = i;
                    CopyVector(xVec, x);
                    return 0;
                }

                var shat = ApplyPreconditioner(preconditioner, s);
                var t = A * shat;
                double normt = t.DotProduct(t);
                if (normt == 0.0)
                {
                    throw new Exception("division by zero");
                }

                omega = t.DotProduct(s) / normt;
                xVec += alpha * phat + omega * shat;
                r = s - omega * t;

                resid = r.L2Norm() / normb;
                if (resid < tol)
                {
                    tol = resid;
                    maxIter = i;
                    CopyVector(xVec, x);
                    return 0;
                }

                if (omega == 0.0)
                {
                    tol = resid;
                    return 3;
                }

                rho2 = rho1;
            }

            tol = resid;
            CopyVector(xVec, x);
            return 1;
        }

        /// <summary>
        /// 无预条件器的 BiCGSTAB。
        /// </summary>
        public static int BiCGSTAB(SparseMatrix A, double[] x, double[] b, ref int maxIter, ref double tol)
        {
            return BiCGSTAB(A, x, b, null, ref maxIter, ref tol);
        }

        /// <summary>
        /// 稀疏上三角矩阵求解。
        /// </summary>
        public static void sparseUpperTriaSolve(SparseMatrix M, double[] x)
        {
            if (M == null)
            {
                throw new Exception("matrix is null");
            }

            if (x == null)
            {
                throw new Exception("vector is null");
            }

            if (M.RowCount != M.ColumnCount)
            {
                throw new Exception("can not use non-quadratic matrix");
            }

            if (M.RowCount != x.Length)
            {
                throw new Exception($"vector size {x.Length} does not match column number {M.ColumnCount}");
            }

            int n = M.RowCount;
            for (int i = n - 1; i >= 0; i--)
            {
                double t = x[i];
                for (int j = i + 1; j < n; j++)
                {
                    double value = M[i, j];
                    if (value != 0.0)
                    {
                        t -= value * x[j];
                    }
                }

                double diag = M[i, i];
                if (diag == 0.0)
                {
                    throw new Exception($"zero diagonal element in row {i}; can not solve linear system");
                }

                x[i] = t / diag;
            }
        }

        /// <summary>
        /// 稀疏下三角矩阵求解。
        /// </summary>
        public static void sparseLowerTriaSolve(SparseMatrix M, double[] x)
        {
            if (M == null)
            {
                throw new Exception("matrix is null");
            }

            if (x == null)
            {
                throw new Exception("vector is null");
            }

            if (M.RowCount != M.ColumnCount)
            {
                throw new Exception("can not use non-quadratic matrix");
            }

            if (M.RowCount != x.Length)
            {
                throw new Exception($"vector size {x.Length} does not match column number {M.ColumnCount}");
            }

            int n = M.RowCount;
            double firstDiag = M[0, 0];
            if (firstDiag == 0.0)
            {
                throw new Exception("zero diagonal element encountered");
            }

            x[0] /= firstDiag;

            for (int i = 1; i < n; i++)
            {
                double t = x[i];
                for (int j = 0; j < i; j++)
                {
                    double value = M[i, j];
                    if (value != 0.0)
                    {
                        t -= value * x[j];
                    }
                }

                double diag = M[i, i];
                if (diag == 0.0)
                {
                    throw new Exception("zero diagonal element encountered");
                }

                x[i] = t / diag;
            }
        }

        /// <summary>
        /// 共轭梯度法求解稀疏线性方程组（不带预条件）。
        /// </summary>
        public static void CGsolve(SparseMatrix a, double[] b, double[] result, double epsilon, ref int kMax)
        {
            if (a == null)
            {
                throw new Exception("matrix is null");
            }

            if (b == null || result == null)
            {
                throw new Exception("vector is null");
            }

            int n = a.RowCount;
            if (a.ColumnCount != n)
            {
                throw new Exception("matrix must be square for CGsolve");
            }

            if (b.Length != n || result.Length != n)
            {
                throw new Exception("vector size does not match matrix size");
            }

            var builder = Vector<double>.Build;
            var r = builder.DenseOfArray(b) - a * builder.DenseOfArray(result);
            var p = r.Clone();
            var w = builder.Dense(n, 0.0);

            double rho = r.DotProduct(r);
            double norm = builder.DenseOfArray(b).DotProduct(builder.DenseOfArray(b));
            double rhoOld = 1.0;

            int k = 0;
            while (rho > epsilon * norm && k < kMax)
            {
                k++;
                if (k > 1)
                {
                    double beta = rho / rhoOld;
                    p = r + beta * p;
                }

                w = a * p;
                double alpha = rho / p.DotProduct(w);
                for (int i = 0; i < n; i++)
                {
                    result[i] += alpha * p[i];
                }

                r -= alpha * w;
                rhoOld = rho;
                rho = r.DotProduct(r);
            }

            if (k == kMax)
            {
                throw new Exception("CGsolve: exceeded maximal number of iterations");
            }

            kMax = k;
        }

        /// <summary>
        /// 预条件器应用（为空则返回输入拷贝）。
        /// </summary>
        private static Vector<double> ApplyPreconditioner(
            Func<Vector<double>, Vector<double>> preconditioner,
            Vector<double> input)
        {
            if (preconditioner == null)
            {
                return input.Clone();
            }

            return preconditioner(input);
        }

        /// <summary>
        /// 将 MathNet 向量拷贝到数组。
        /// </summary>
        private static void CopyVector(Vector<double> source, double[] target)
        {
            for (int i = 0; i < target.Length; i++)
            {
                target[i] = source[i];
            }
        }
    }
}
