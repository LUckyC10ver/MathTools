using System;

namespace MathTools.Core
{
    public partial class Functions
    {
        /// <summary>
        /// LU 分解（带主元选取）。
        /// </summary>
        public static void LUDecompose(double[][] a, int[] indx, out double d)
        {
            if (a == null || a.Length == 0)
            {
                throw new Exception("matrix is empty");
            }

            int n = a.Length;
            if (a[0].Length != n)
            {
                throw new Exception("matrix must be square for LU decomposition");
            }

            if (indx == null || indx.Length != n)
            {
                throw new Exception("permutation vector has wrong size");
            }

            const double tiny = 1e-30;
            var vv = new double[n];
            d = 1.0;

            for (int i = 0; i < n; i++)
            {
                double big = 0.0;
                for (int j = 0; j < n; j++)
                {
                    double temp = Math.Abs(a[i][j]);
                    if (temp > big)
                    {
                        big = temp;
                    }
                }

                if (big == 0.0)
                {
                    throw new Exception("Singular matrix");
                }

                vv[i] = 1.0 / big;
            }

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < j; i++)
                {
                    double sum = a[i][j];
                    for (int k = 0; k < i; k++)
                    {
                        sum -= a[i][k] * a[k][j];
                    }

                    a[i][j] = sum;
                }

                double big = 0.0;
                int imax = 0;
                for (int i = j; i < n; i++)
                {
                    double sum = a[i][j];
                    for (int k = 0; k < j; k++)
                    {
                        sum -= a[i][k] * a[k][j];
                    }

                    a[i][j] = sum;
                    double dum = vv[i] * Math.Abs(sum);
                    if (Math.Abs(dum) >= Math.Abs(big))
                    {
                        big = dum;
                        imax = i;
                    }
                }

                if (j != imax)
                {
                    for (int k = 0; k < n; k++)
                    {
                        double dum = a[imax][k];
                        a[imax][k] = a[j][k];
                        a[j][k] = dum;
                    }

                    d = -d;
                    vv[imax] = vv[j];
                }

                indx[j] = imax;
                if (a[j][j] == 0.0)
                {
                    a[j][j] = tiny;
                }

                if (j != n - 1)
                {
                    double dum = 1.0 / a[j][j];
                    for (int i = j + 1; i < n; i++)
                    {
                        a[i][j] *= dum;
                    }
                }
            }
        }

        /// <summary>
        /// LU 回代求解。
        /// </summary>
        public static void LUBacksubstitute(double[][] a, int[] indx, double[] b)
        {
            if (a == null || a.Length == 0)
            {
                throw new Exception("matrix is empty");
            }

            int n = a.Length;
            if (a[0].Length != n)
            {
                throw new Exception("matrix must be square for LU decomposition");
            }

            if (b == null || b.Length != n)
            {
                throw new Exception("vector size does not match matrix size");
            }

            if (indx == null || indx.Length != n)
            {
                throw new Exception("permutation vector has wrong size");
            }

            int ii = -1;
            for (int i = 0; i < n; i++)
            {
                int ip = indx[i];
                double sum = b[ip];
                b[ip] = b[i];
                if (ii >= 0)
                {
                    for (int j = ii; j <= i - 1; j++)
                    {
                        sum -= a[i][j] * b[j];
                    }
                }
                else if (Math.Abs(sum) > 0.0)
                {
                    ii = i;
                }

                b[i] = sum;
            }

            for (int i = n - 1; i >= 0; i--)
            {
                double sum = b[i];
                for (int j = i + 1; j < n; j++)
                {
                    sum -= a[i][j] * b[j];
                }

                b[i] = sum / a[i][i];
            }
        }
    }
}
