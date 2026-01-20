using System;

namespace MathTools.Core
{
    public partial class Functions
    {
        public static void linprog(
            out double[] x,
            out int[] actSet,
            double[] f,
            double[][] A,
            double[] b,
            double[][] Aeq,
            double[] beq,
            double EPS = 1e-16)
        {
            if (f == null)
            {
                throw new Exception("objective is null");
            }

            int nvars = f.Length;
            int nieq = b?.Length ?? 0;
            int neq = beq?.Length ?? 0;
            if (nieq > 0 && A == null)
            {
                throw new Exception("inequality matrix is null");
            }

            if (neq > 0 && Aeq == null)
            {
                throw new Exception("equality matrix is null");
            }
            int rows = nieq + neq + 2;
            int cols = 2 * nvars + 1;
            double[][] tab = CreateMatrix(rows, cols, 0.0);

            for (int i = 0; i < nvars; i++)
            {
                tab[0][i + 1] = -f[i];
                tab[0][i + nvars + 1] = f[i];
            }

            int icnt = 0;
            int[] index = new int[nieq];
            int k = 1;

            for (int i = 0; i < nieq; i++)
            {
                if (b[i] >= 0.0)
                {
                    tab[k][0] = b[i];
                    for (int j = 0; j < nvars; j++)
                    {
                        double value = A[i][j];
                        tab[k][j + 1] = -value;
                        tab[k][j + nvars + 1] = value;
                    }
                    k++;
                }
                else
                {
                    index[icnt++] = i;
                }
            }

            for (int i = 0; i < icnt; i++)
            {
                tab[k][0] = -b[index[i]];
                for (int j = 0; j < nvars; j++)
                {
                    double value = A[index[i]][j];
                    tab[k][j + nvars + 1] = -value;
                    tab[k][j + 1] = value;
                }
                k++;
            }

            for (int i = 0; i < neq; i++)
            {
                double sgn = beq[i] >= 0.0 ? 1.0 : -1.0;
                tab[nieq + i + 1][0] = sgn * beq[i];
                for (int j = 0; j < nvars; j++)
                {
                    double value = sgn * Aeq[i][j];
                    tab[nieq + i + 1][j + 1] = -value;
                    tab[nieq + i + 1][j + nvars + 1] = value;
                }
            }

            int retVal = 0;
            simplex(tab, nieq - icnt, icnt, neq, ref retVal, out int[] izrov, out int[] iposv, EPS);

            if (retVal == -1)
            {
                throw new Exception("no feasible point");
            }

            if (retVal == 1)
            {
                throw new Exception("objective is unbounded");
            }

            x = new double[nvars];
            for (int j = 1; j <= iposv.Length; j++)
            {
                int idx = iposv[j - 1] - 1;
                if (idx < nvars)
                {
                    x[idx] += tab[j][0];
                }
                else if (idx < 2 * nvars)
                {
                    x[idx - nvars] -= tab[j][0];
                }
            }

            actSet = new int[0];
            if (neq > 0)
            {
                actSet = new int[neq];
                for (int i = 0; i < neq; i++)
                {
                    actSet[i] = i;
                }
            }

            if (nieq > 0)
            {
                double[] resIeq = matdotvecplusvec(A, x, b, -1.0);
                int count = actSet.Length;
                Array.Resize(ref actSet, count + resIeq.Length);
                for (int i = 0; i < resIeq.Length; i++)
                {
                    if (Math.Abs(resIeq[i]) < EPS)
                    {
                        actSet[count++] = i + neq;
                    }
                }

                Array.Resize(ref actSet, count);
            }
        }

        public static void simplex(
            double[][] a,
            int m1,
            int m2,
            int m3,
            ref int icase,
            out int[] izrov,
            out int[] iposv,
            double EPS = 1e-16)
        {
            int m = a.Length - 2;
            int n = a[0].Length - 1;
            if (m != (m1 + m2 + m3))
            {
                throw new Exception("Bad input constraint counts in simplex");
            }

            int[] l1 = new int[n + 1];
            int[] l2 = new int[m];
            int[] l3 = new int[m];

            izrov = new int[n];
            iposv = new int[m];

            int nl1 = n;
            for (int k = 1; k <= n; k++)
            {
                l1[k - 1] = izrov[k - 1] = k;
            }

            int nl2 = m;
            for (int i = 1; i <= m; i++)
            {
                if (a[i][0] < 0.0)
                {
                    throw new Exception("Bad input tableau in simplex");
                }

                l2[i - 1] = i;
                iposv[i - 1] = n + i;
            }

            for (int i = 1; i <= m2; i++)
            {
                l3[i - 1] = 1;
            }

            int ir = 0;
            if (m2 + m3 > 0)
            {
                ir = 1;
                for (int k = 1; k <= (n + 1); k++)
                {
                    double q1 = 0.0;
                    for (int i = m1 + 1; i <= m; i++)
                    {
                        q1 += a[i][k - 1];
                    }
                    a[m + 1][k - 1] = -q1;
                }

                do
                {
                    simp1(a, m + 1, l1, nl1, 0, out int kp, out double bmax);

                    if (bmax <= EPS && a[m + 1][0] < -EPS)
                    {
                        icase = -1;
                        return;
                    }

                    if (bmax <= EPS && a[m + 1][0] <= EPS)
                    {
                        int m12 = m1 + m2 + 1;
                        if (m12 <= m)
                        {
                            for (int ip = m12; ip <= m; ip++)
                            {
                                if (iposv[ip - 1] == (ip + n))
                                {
                                    simp1(a, ip, l1, nl1, 1, out kp, out bmax);
                                    if (bmax > EPS)
                                    {
                                        goto one;
                                    }
                                }
                            }
                        }

                        ir = 0;
                        m12--;
                        if (m1 + 1 <= m12)
                        {
                            for (int i = m1 + 1; i <= m12; i++)
                            {
                                if (l3[i - m1 - 1] == 1)
                                {
                                    for (int k = 1; k <= n + 1; k++)
                                    {
                                        a[i][k - 1] = -a[i][k - 1];
                                    }
                                }
                            }
                        }
                        break;
                    }

                    simp2(a, n, l2, nl2, out int ip2, kp, out double q1b, EPS);
                    if (ip2 == 0)
                    {
                        icase = -1;
                        return;
                    }

                one:
                    simp3(a, m + 1, n, ip2, kp);

                    if (iposv[ip2 - 1] >= (n + m1 + m2 + 1))
                    {
                        int kk = 0;
                        for (int k = 1; k <= nl1; k++)
                        {
                            if (l1[k - 1] == kp)
                            {
                                kk = k;
                                break;
                            }
                        }

                        nl1--;
                        for (int is = kk; is <= nl1; is++)
                        {
                            l1[is - 1] = l1[is];
                        }

                        a[m + 1][kp] += 1.0;
                        for (int i = 1; i <= m + 2; i++)
                        {
                            a[i - 1][kp] = -a[i - 1][kp];
                        }
                    }
                    else if (iposv[ip2 - 1] >= (n + m1 + 1))
                    {
                        int kh = iposv[ip2 - 1] - m1 - n;
                        if (l3[kh - 1] == 1)
                        {
                            l3[kh - 1] = 0;
                            a[m + 1][kp] += 1.0;
                            for (int i = 1; i <= m + 2; i++)
                            {
                                a[i - 1][kp] = -a[i - 1][kp];
                            }
                        }
                    }

                    int is2 = izrov[kp - 1];
                    izrov[kp - 1] = iposv[ip2 - 1];
                    iposv[ip2 - 1] = is2;
                } while (ir != 0);
            }

            for (;;)
            {
                simp1(a, 0, l1, nl1, 0, out int kp, out double bmax);
                if (bmax < EPS)
                {
                    icase = 0;
                    return;
                }

                simp2(a, n, l2, nl2, out int ip, kp, out double q1, EPS);
                if (ip == 0)
                {
                    icase = 1;
                    return;
                }

                simp3(a, m, n, ip, kp);

                int is = izrov[kp - 1];
                izrov[kp - 1] = iposv[ip - 1];
                iposv[ip - 1] = is;
            }
        }

        private static void simp1(double[][] a, int mm, int[] ll, int nll, int iabf, out int kp, out double bmax)
        {
            kp = ll[0];
            bmax = a[mm][kp];
            for (int k = 2; k <= nll; k++)
            {
                double test;
                if (iabf == 0)
                {
                    test = a[mm][ll[k - 1]] - bmax;
                }
                else
                {
                    test = Math.Abs(a[mm][ll[k - 1]]) - Math.Abs(bmax);
                }

                if (test > 0.0)
                {
                    bmax = a[mm][ll[k - 1]];
                    kp = ll[k - 1];
                }
            }
        }

        private static void simp2(double[][] a, int n, int[] l2, int nl2, out int ip, int kp, out double q1, double EPS)
        {
            ip = 0;
            q1 = 0.0;
            for (int i = 1; i <= nl2; i++)
            {
                if (a[l2[i - 1]][kp] < -EPS)
                {
                    q1 = -a[l2[i - 1]][0] / a[l2[i - 1]][kp];
                    ip = l2[i - 1];
                    for (int ii = i + 1; ii <= nl2; ii++)
                    {
                        int row = l2[ii - 1];
                        if (a[row][kp] < -EPS)
                        {
                            double q = -a[row][0] / a[row][kp];
                            if (q < q1)
                            {
                                ip = row;
                                q1 = q;
                            }
                            else if (q == q1)
                            {
                                double qp = 0.0;
                                double q0 = 0.0;
                                for (int k = 1; k <= n; k++)
                                {
                                    qp = -a[ip][k] / a[ip][kp];
                                    q0 = -a[row][k] / a[row][kp];
                                    if (q0 != qp)
                                    {
                                        break;
                                    }
                                }
                                if (q0 < qp)
                                {
                                    ip = row;
                                }
                            }
                        }
                    }
                    break;
                }
            }
        }

        private static void simp3(double[][] a, int i1, int k1, int ip, int kp)
        {
            double piv = 1.0 / a[ip][kp];
            for (int ii = 1; ii <= i1 + 1; ii++)
            {
                if (ii - 1 != ip)
                {
                    a[ii - 1][kp] *= piv;
                    for (int kk = 1; kk <= k1 + 1; kk++)
                    {
                        if (kk - 1 != kp)
                        {
                            a[ii - 1][kk - 1] -= a[ip][kk - 1] * a[ii - 1][kp];
                        }
                    }
                }
            }

            for (int kk = 1; kk <= k1 + 1; kk++)
            {
                if (kk - 1 != kp)
                {
                    a[ip][kk - 1] *= -piv;
                }
            }

            a[ip][kp] = piv;
        }

        private static double[][] CreateMatrix(int rows, int cols, double value)
        {
            var matrix = new double[rows][];
            for (int i = 0; i < rows; i++)
            {
                matrix[i] = new double[cols];
                if (value != 0.0)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        matrix[i][j] = value;
                    }
                }
            }

            return matrix;
        }
    }
}
