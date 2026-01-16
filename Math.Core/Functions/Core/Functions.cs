using System;
using System.Collections.Generic;
using System.Linq;
using MathTools.Core.Extensions;
using MathNet.Numerics.LinearAlgebra;

namespace MathTools.Core
{
    public partial class Functions
    {
        /// <summary>
        /// 求中位数
        /// </summary>
        /// <param name="input"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static double median(double[] input)
        {
            if (input == null || input.Length == 0)
            {
                throw new Exception("empty or irregular vector");
            }

            var sorted = (double[])input.Clone();
            Array.Sort(sorted);

            int count = sorted.Length;
            if (count % 2 == 0)
            {
                int upper = count / 2;
                int lower = upper - 1;
                return (sorted[lower] + sorted[upper]) / 2.0;
            }

            return sorted[count / 2];
        }

        /// <summary>
        /// 一维数组绝对值最大值
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static double max_abs(double[] a)
        {
            if (a == null || a.Length == 0)
            {
                return 0;
            }

            return a.ToVector().AbsoluteMaximum();
        }

        /// <summary>
        /// 一维数组范数
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static double norm(double[] a)
        {
            return a.ToVector().L2Norm();
        }

        /// <summary>
        /// 霍纳法则，求解多项式
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double horner(double[] a, double x)
        {
            if (a == null || a.Length == 0)
            {
                return 0;
            }

            double result = a[0];
            for (int i = 1; i < a.Length; i++)
            {
                result = result * x + a[i];
            }

            return result;
        }

        /// <summary>
        /// 计算多项式在给定点 x 处的导数值
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double hornerDerivative(double[] a, double x)
        {
            if (a == null || a.Length == 0)
            {
                return 0;
            }

            int grad = a.Length - 1;
            if (grad == 0)
            {
                return 0;
            }

            double result = grad * a[0];
            for (int i = 1; i < grad; i++)
            {
                result = result * x + (grad - i) * a[i];
            }

            return result;
        }

        /// <summary>
        /// 计算当前数组的百分位数组
        /// </summary>
        /// <param name="input"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static double[] quantile(double[] input, double[] p)
        {
            if (input == null || input.Length == 0)
            {
                throw new Exception("empty or irregular vector");
            }

            var sorted = (double[])input.Clone();
            Array.Sort(sorted);

            int noPoints = sorted.Length;
            double[] result = new double[p.Length];
            for (int i = 0; i < p.Length; i++)
            {
                if (p[i] < 0 || p[i] > 1)
                {
                    throw new Exception("entry of p is not in the interval [0 1]");
                }

                double nonIntIndex = p[i] * (noPoints - 1);
                int intPart = (int)Math.Floor(nonIntIndex);
                double fractPart = nonIntIndex - intPart;
                int lowerIndex = intPart;
                int upperIndex = lowerIndex + 1;
                if ((lowerIndex == noPoints - 1) || (fractPart < 1e-15))
                {
                    result[i] = sorted[lowerIndex];
                }
                else
                {
                    result[i] = sorted[lowerIndex] * (1.0 - fractPart) + sorted[upperIndex] * fractPart;
                }
            }

            return result;
        }

        /// <summary>
        /// 一维数组点乘
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static double dot(double[] a, double[] b)
        {
            if (a == null || b == null || a.Length != b.Length)
            {
                throw new Exception("vector sizes differ in dot<VECdoubleOR>( VECdoubleOR&, VECdoubleOR&)");
            }

            return a.ToVector().DotProduct(b.ToVector());
        }

        /// <summary>
        /// 一维数组字典序比较
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="tol"></param>
        /// <returns></returns>
        public static bool less_lexicographical(double[] a, double[] b, double tol = 0)
        {
            int n = Math.Min(a.Length, b.Length);

            for (int i = 0; i < n; i++)
            {
                if (a[i] + tol < b[i])
                {
                    return true;
                }

                if (b[i] + tol < a[i])
                {
                    return false;
                }
            }

            return a.Length < b.Length;
        }

        /// <summary>
        /// 希尔排序算法
        /// </summary>
        /// <param name="xValues"></param>
        /// <param name="yValues"></param>
        public static void sortPoints(ref double[] xValues, ref double[] yValues)
        {
            if (xValues == null || yValues == null || xValues.Length != yValues.Length)
            {
                throw new Exception("vector sizes differ in sortPoints");
            }

            int n = xValues.Length;
            if (n < 2)
            {
                return;
            }

            int inc = 1;
            while (inc < n)
            {
                inc = inc * 3 + 1;
            }

            while (inc > 1)
            {
                inc /= 3;
                for (int i = inc; i < n; i++)
                {
                    double v = xValues[i];
                    double y = yValues[i];
                    int j = i;
                    while (j >= inc && xValues[j - inc] > v)
                    {
                        xValues[j] = xValues[j - inc];
                        yValues[j] = yValues[j - inc];
                        j -= inc;
                    }

                    xValues[j] = v;
                    yValues[j] = y;
                }
            }
        }

        /// <summary>
        /// 数组标准差
        /// </summary>
        /// <param name="aVector"></param>
        /// <returns></returns>
        public static double vecStddev(double[] aVector)
        {
            if (aVector == null || aVector.Length == 0)
            {
                return 0;
            }

            double sumX = 0;
            double sumX2 = 0;

            for (var i = 0; i < aVector.Length; i++)
            {
                sumX += aVector[i];
                sumX2 += sqrInternal(aVector[i]);
            }

            if (aVector.Length > 1)
            {
                double help = sumX2 - sqrInternal(sumX) / aVector.Length;
                help /= aVector.Length - 1;
                if (help < 0)
                {
                    help = 0;
                }

                return Math.Sqrt(help);
            }

            return 0;
        }

        /// <summary>
        /// 返回数组最大值索引
        /// </summary>
        /// <param name="vec"></param>
        /// <param name="ind"></param>
        /// <returns></returns>
        public static double getMaximum(double[] vec, out int ind)
        {
            if (vec == null || vec.Length == 0)
            {
                ind = -1;
                return 0;
            }

            double temp = Math.Abs(vec[0]);
            int maxi = 0;

            for (int i = 1; i < vec.Length; i++)
            {
                double abs = Math.Abs(vec[i]);
                if (abs > temp)
                {
                    temp = abs;
                    maxi = i;
                }
            }

            ind = maxi;
            return temp;
        }

        /// <summary>
        /// 对一维数组进行重新排列或打乱顺序
        /// </summary>
        /// <param name="v"></param>
        /// <param name="perm"></param>
        /// <exception cref="Exception"></exception>
        public static void permuteVector(ref double[] v, int[] perm)
        {
            int n = v.Length;
            int n1 = perm.Length;

            if (n != n1)
            {
                throw new Exception("size(v) != size(perm)");
            }

            double[] copy = (double[])v.Clone();

            for (int i = 0; i < n; i++)
            {
                int pi = perm[i];
                if ((pi < 0) || (pi >= n))
                {
                    throw new Exception("badly chosen permutation");
                }

                v[pi] = copy[i];
            }
        }

        /// <summary>
        /// 初始化一个一维线性空间
        /// </summary>
        /// <param name="vec"></param>
        /// <param name="x1"></param>
        /// <param name="x2"></param>
        /// <param name="count"></param>
        public static double[] linspace(double x1, double x2, int count)
        {
            var vec = new double[count];
            if (count > 1)
            {
                var deltaX = (x2 - x1) / (count - 1);
                for (int i = 0; i < count; i++)
                {
                    vec[i] = x1 + i * deltaX;
                }
            }

            return vec;
        }

        /// <summary>
        /// vec重新排序，返回排序后的索引
        /// </summary>
        /// <param name="vec"></param>
        public static int[] sort(double[] vec)
        {
            int n = vec.Length;
            int[] index = new int[n];
            for (int i = 0; i < n; i++)
            {
                index[i] = i;
            }

            var values = (double[])vec.Clone();
            Array.Sort(values, index);

            return index;
        }

        /// <summary>
        /// 两个整数列表的交集查找
        /// </summary>
        /// <param name="vec1"></param>
        /// <param name="vec2"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static int[] intersect(int[] vec1, int[] vec2)
        {
            var v1 = (int[])vec1.Clone();
            var v2 = (int[])vec2.Clone();
            Array.Sort(v1);
            Array.Sort(v2);

            int iv1 = 0;
            int iv2 = 0;
            var result = new List<int>();
            while ((iv1 < v1.Length) && (iv2 < v2.Length))
            {
                if (v1[iv1] == v2[iv2])
                {
                    result.Add(v1[iv1]);
                    iv1++;
                    iv2++;
                }
                else if (v1[iv1] < v2[iv2])
                {
                    iv1++;
                }
                else if (v2[iv2] < v1[iv1])
                {
                    iv2++;
                }
                else
                {
                    throw new Exception("wrong definition of operator<");
                }
            }

            return result.ToArray();
        }

        /// <summary>
        /// 对一个矩阵 a 进行秩1更新，即通过向量 x 和权重 weight 来调整矩阵 a 的值
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <param name="weight"></param>
        /// <exception cref="Exception"></exception>
        public static void rank1Update(ref double[][] a, double[] x, double weight)
        {
            int rows = a.Length;
            int cols = a[0].Length;
            if (x.Length != rows || x.Length != cols)
            {
                throw new Exception("incompatible matrix/vector sizes");
            }

            var v = x.ToVector();
            var updated = a.ToMatrix() + v.OuterProduct(v) * weight;
            ReplaceMatrix(ref a, updated);
        }

        /// <summary>
        /// 计算矩阵的迹
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static double trace(double[][] a)
        {
            return a.ToMatrix().Trace();
        }

        /// <summary>
        /// 矩阵点乘
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static double[][] matdotmat(double[][] left, double[][] right)
        {
            return left.ToMatrix().Multiply(right.ToMatrix()).ToJaggedArray();
        }

        /// <summary>
        /// 矩阵转置点乘
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] transmatdotmat(double[][] left, double[][] right)
        {
            return left.ToMatrix().Transpose().Multiply(right.ToMatrix()).ToJaggedArray();
        }

        /// <summary>
        /// 矩阵乘向量配合加法
        /// </summary>
        /// <param name="A"></param>
        /// <param name="x"></param>
        /// <param name="b"></param>
        /// <param name="s"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static double[] matdotvecplusvec(double[][] A, double[] x, double[] b, double s = 1.0f)
        {
            var result = A.ToMatrix().Multiply(x.ToVector()) + b.ToVector() * s;
            return result.ToArray();
        }

        /// <summary>
        /// 将矩阵 H 投影到矩阵 Z 上
        /// </summary>
        /// <param name="H"></param>
        /// <param name="Z"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static double[][] projectMatrix(double[][] H, double[][] Z)
        {
            var h = H.ToMatrix();
            var z = Z.ToMatrix();
            if (h.ColumnCount != z.RowCount || h.RowCount != h.ColumnCount)
            {
                throw new Exception("incompatible matrix sizes");
            }

            return z.TransposeThisAndMultiply(h).Multiply(z).ToJaggedArray();
        }

        /// <summary>
        /// 计算矩阵最长边
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static int length_m(double[][] m)
        {
            return Math.Max(m.Length, m[0].Length);
        }

        /// <summary>
        /// 矩阵所有元素写1
        /// </summary>
        /// <param name="m1"></param>
        public static void ones_matrix(ref double[][] m1)
        {
            ReplaceMatrix(ref m1, Matrix<double>.Build.Dense(m1.Length, m1[0].Length, 1.0));
        }

        /// <summary>
        /// 矩阵所有元素乘val
        /// </summary>
        /// <param name="m"></param>
        /// <param name="val"></param>
        /// <returns></returns>
        public static void mul_m_s(ref double[][] m, double val)
        {
            ReplaceMatrix(ref m, m.ToMatrix() * val);
        }

        /// <summary>
        /// 矩阵所有元素求绝对值
        /// </summary>
        /// <param name="matrix"></param>
        public static void abs_m(ref double[][] matrix)
        {
            ReplaceMatrix(ref matrix, matrix.ToMatrix().Map(Math.Abs));
        }

        /// <summary>
        /// 矩阵相加
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] add_m_mm(double[][] left, double[][] right)
        {
            return left.ToMatrix().Add(right.ToMatrix()).ToJaggedArray();
        }

        /// <summary>
        /// 矩阵所有元素写0
        /// </summary>
        /// <param name="m1"></param>
        public static void zero_matrix(ref double[][] m1)
        {
            ReplaceMatrix(ref m1, Matrix<double>.Build.Dense(m1.Length, m1[0].Length, 0.0));
        }

        /// <summary>
        /// 构造单位矩阵
        /// </summary>
        /// <param name="m1"></param>
        /// <returns></returns>
        public static void unit_matrix(ref double[][] m1)
        {
            ReplaceMatrix(ref m1, Matrix<double>.Build.DenseIdentity(m1.Length, m1[0].Length));
        }

        /// <summary>
        /// 矩阵数乘输出新的矩阵
        /// </summary>
        /// <param name="src"></param>
        /// <param name="val"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] mul_m_ms(double[][] src, double val)
        {
            return (src.ToMatrix() * val).ToJaggedArray();
        }

        /// <summary>
        /// 矩阵除法，输出新的矩阵
        /// </summary>
        /// <param name="val"></param>
        /// <param name="src"></param>
        /// <param name="EPS"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] div_m_sm(double val, double[][] src, double EPS = 1.0e-20)
        {
            if (Math.Abs(val) < EPS)
            {
                throw new Exception("division by zero");
            }

            return mul_m_ms(src, 1 / val);
        }

        /// <summary>
        /// 矩阵求绝对值，输出新的矩阵
        /// </summary>
        /// <param name="src"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] abs_m_m(double[][] src)
        {
            return src.ToMatrix().Map(Math.Abs).ToJaggedArray();
        }

        /// <summary>
        /// 矩阵求符号位，输出新的矩阵
        /// </summary>
        /// <param name="src"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] sign_m_m(double[][] src)
        {
            return src.ToMatrix().Map(value => value > 0 ? 1.0 : -1.0).ToJaggedArray();
        }

        /// <summary>
        /// 矩阵乘法，输出新的矩阵
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] mul_m_mm(double[][] left, double[][] right)
        {
            return left.ToMatrix().PointwiseMultiply(right.ToMatrix()).ToJaggedArray();
        }

        /// <summary>
        /// 矩阵加法，输出新的矩阵
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] sub_m_mm(double[][] left, double[][] right)
        {
            return left.ToMatrix().Subtract(right.ToMatrix()).ToJaggedArray();
        }

        /// <summary>
        /// 矩阵转置，输出新的矩阵
        /// </summary>
        /// <param name="src"></param>
        public static double[][] transpose(double[][] src)
        {
            return src.ToMatrix().Transpose().ToJaggedArray();
        }

        /// <summary>
        /// 矩阵除法
        /// </summary>
        /// <param name="m"></param>
        /// <param name="val"></param>
        /// <param name="EPS"></param>
        /// <exception cref="Exception"></exception>
        public static void div_m_s(ref double[][] m, double val, double EPS = 1.0e-20)
        {
            if (Math.Abs(val) < EPS)
            {
                throw new Exception("division by zero");
            }

            ReplaceMatrix(ref m, m.ToMatrix() / val);
        }

        /// <summary>
        /// 矩阵有限差分法（差输出单元，源代码输出字节流）
        /// </summary>
        /// <param name="finite_diff_deriv"></param>
        /// <param name="analytic_deriv"></param>
        /// <param name="evalstr2"></param>
        /// <param name="outx"></param>
        /// <exception cref="Exception"></exception>
        public static void graderr_m(double[][] finite_diff_deriv, double[][] analytic_deriv, string evalstr2, byte[] outx)
        {
            int rows = finite_diff_deriv.Length;
            int cols = finite_diff_deriv[0].Length;
            if (analytic_deriv.Length != rows || analytic_deriv[0].Length != cols)
            {
                throw new Exception("finite_diff_deriv and analytic_deriv are incompatible");
            }

            var help_1 = sub_m_mm(finite_diff_deriv, analytic_deriv);
            abs_m(ref help_1);
            double err = max_val(help_1, 1.0e+99);

            double nor = normFrobenius(analytic_deriv);
            if (err > 1e-6 * nor + 1e-5)
            {
                return;
            }
        }

        /// <summary>
        /// 从大矩阵中提取一个小矩阵
        /// </summary>
        /// <param name="m"></param>
        /// <param name="r_low"></param>
        /// <param name="r_high"></param>
        /// <param name="c_low"></param>
        /// <param name="c_high"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] submat(double[][] m, int r_low, int r_high, int c_low, int c_high)
        {
            int ex_rows = r_high - r_low + 1;
            int ex_cols = c_high - c_low + 1;
            return m.ToMatrix().SubMatrix(r_low, ex_rows, c_low, ex_cols).ToJaggedArray();
        }

        /// <summary>
        /// 输出该矩阵最小值
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="MAX"></param>
        /// <returns></returns>
        public static double min_val(double[][] matrix, double MAX)
        {
            var min = matrix.ToMatrix().Enumerate().Min();
            return Math.Min(min, MAX);
        }

        /// <summary>
        /// 输出无穷范数或者矩阵行最大值
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static double normInfinity(double[][] m)
        {
            return m.ToMatrix().LInfinityNorm();
        }

        /// <summary>
        /// 计算矩阵的Frobenius范数
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static double normFrobenius(double[][] m)
        {
            return m.ToMatrix().FrobeniusNorm();
        }

        /// <summary>
        /// 求矩阵最大值
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="MAX"></param>
        /// <returns></returns>
        public static double max_val(double[][] matrix, double MAX)
        {
            var max = matrix.ToMatrix().Enumerate().Max();
            return Math.Max(max, -MAX);
        }

        /// <summary>
        /// 矩阵加法
        /// </summary>
        /// <param name="target"></param>
        /// <param name="src"></param>
        /// <exception cref="Exception"></exception>
        public static void add_m_m(ref double[][] target, double[][] src)
        {
            ReplaceMatrix(ref target, target.ToMatrix().Add(src.ToMatrix()));
        }

        /// <summary>
        /// 矩阵减法
        /// </summary>
        /// <param name="target"></param>
        /// <param name="src"></param>
        /// <exception cref="Exception"></exception>
        public static void sub_m_m(ref double[][] target, double[][] src)
        {
            ReplaceMatrix(ref target, target.ToMatrix().Subtract(src.ToMatrix()));
        }

        /// <summary>
        /// 矩阵元素乘法
        /// </summary>
        /// <param name="target"></param>
        /// <param name="src"></param>
        /// <exception cref="Exception"></exception>
        public static void mul_m_m(ref double[][] target, double[][] src)
        {
            ReplaceMatrix(ref target, target.ToMatrix().PointwiseMultiply(src.ToMatrix()));
        }

        /// <summary>
        /// 矩阵加常数
        /// </summary>
        /// <param name="m"></param>
        /// <param name="val"></param>
        public static void add_m_s(ref double[][] m, double val)
        {
            ReplaceMatrix(ref m, m.ToMatrix().Map(value => value + val));
        }

        /// <summary>
        /// 矩阵元素除法
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <param name="EPS"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] div_m_mm(double[][] left, double[][] right, double EPS)
        {
            var rightMatrix = right.ToMatrix();
            var minAbs = rightMatrix.Enumerate().Select(Math.Abs).DefaultIfEmpty(0.0).Min();
            if (sqrInternal(minAbs) < sqrInternal(EPS))
            {
                throw new Exception("division by zero");
            }

            return left.ToMatrix().PointwiseDivide(rightMatrix).ToJaggedArray();
        }

        /// <summary>
        /// 将一个矩阵插入另一个矩阵
        /// </summary>
        /// <param name="t"></param>
        /// <param name="s"></param>
        /// <param name="row"></param>
        /// <param name="col"></param>
        /// <exception cref="Exception"></exception>
        public static void set_m_ss(ref double[][] t, double[][] s, int row, int col)
        {
            int t_rows = t.Length;
            int t_cols = t[0].Length;
            if ((row < 0) || (col < 0) || (row >= t_rows) || (col >= t_cols))
            {
                throw new Exception("illegal row or col-values!");
            }

            int s_rows = s.Length;
            int s_cols = s[0].Length;
            if ((row + s_rows > t_rows) || (col + s_cols > t_cols))
            {
                throw new Exception("capacity of target matrix exceeded!");
            }

            for (int i = 0; i < s_rows; i++)
            {
                for (int j = 0; j < s_cols; j++)
                {
                    t[i + row][j + col] = s[i][j];
                }
            }
        }

        /// <summary>
        /// 矩阵元素加法，输出新的矩阵
        /// </summary>
        /// <param name="src"></param>
        /// <param name="val"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] add_m_ms(double[][] src, double val)
        {
            return src.ToMatrix().Map(value => value + val).ToJaggedArray();
        }

        /// <summary>
        /// 矩阵对角线组合
        /// </summary>
        /// <param name="t"></param>
        /// <param name="m1"></param>
        /// <param name="m1rl"></param>
        /// <param name="m1rh"></param>
        /// <param name="m2"></param>
        /// <param name="m2rl"></param>
        /// <param name="m2rh"></param>
        /// <exception cref="Exception"></exception>
        public static void connect_rows(ref double[][] t, double[][] m1, int m1rl, int m1rh, double[][] m2, int m2rl, int m2rh)
        {
            int m1_cols = m1[0].Length;
            int m2_cols = m2[0].Length;

            if (m1_cols != m2_cols)
            {
                throw new Exception("Can't connect matrices, #columns are different!");
            }

            int m1r = m1rl > m1rh ? 0 : m1rh - m1rl + 1;
            int m2r = m2rl > m2rh ? 0 : m2rh - m2rl + 1;
            if ((t.Length != m1r + m2r) || (t[0].Length != m1_cols))
            {
                throw new Exception("target has wrong Dimensions!");
            }

            for (int i = 0; i < m1r; i++)
            {
                for (int j = 0; j < m1_cols; j++)
                {
                    t[i][j] = m1[m1rl + i][j];
                }
            }

            for (int i = 0; i < m2r; i++)
            {
                for (int j = 0; j < m2_cols; j++)
                {
                    t[i + m1r][j] = m2[m2rl + i][j];
                }
            }
        }

        /// <summary>
        /// 求矩阵列最小值
        /// </summary>
        /// <param name="rowvec"></param>
        /// <param name="res"></param>
        /// <param name="src"></param>
        /// <param name="MAX"></param>
        /// <exception cref="Exception"></exception>
        public static void mincols(ref double[][] rowvec, ref int[] res, double[][] src, double MAX)
        {
            int rows = src.Length;
            int cols = src[0].Length;
            if ((rowvec.Length != 1) || (rowvec[0].Length != cols))
            {
                throw new Exception("incorrect dimension for result vector");
            }

            if (res.Length != cols)
            {
                throw new Exception("incorrect dimension for index result vector");
            }

            for (int j = 0; j < cols; j++)
            {
                double min = MAX;
                int min_index = 0;
                for (int i = 0; i < rows; i++)
                {
                    double el = src[i][j];
                    if (el < min)
                    {
                        min = el;
                        min_index = i;
                    }
                }

                rowvec[0][j] = min;
                res[j] = min_index;
            }
        }

        /// <summary>
        /// 求矩阵行最小值
        /// </summary>
        /// <param name="colvec"></param>
        /// <param name="res"></param>
        /// <param name="src"></param>
        /// <param name="MAX"></param>
        /// <exception cref="Exception"></exception>
        public static void minrows(ref double[][] colvec, ref int[] res, double[][] src, double MAX)
        {
            int rows = src.Length;
            int cols = src[0].Length;

            if ((colvec.Length != rows) || (colvec[0].Length != 1))
            {
                throw new Exception("MINROWS: incorrect dimension for result vector");
            }

            if (res.Length != rows)
            {
                throw new Exception("MINROWS: incorrect dimension for index result vector");
            }

            for (int i = 0; i < rows; i++)
            {
                double min = MAX;
                int min_index = 0;

                for (int j = 0; j < cols; j++)
                {
                    double el = src[i][j];
                    if (el < min)
                    {
                        min = el;
                        min_index = j;
                    }
                }

                colvec[i][0] = min;
                res[i] = min_index;
            }
        }

        /// <summary>
        /// 求最小值矩阵
        /// </summary>
        /// <param name="rowvec"></param>
        /// <param name="res"></param>
        /// <param name="src"></param>
        public static void min_mat(ref double[][] rowvec, ref int[] res, double[][] src)
        {
            if ((src.Length > 1) && (src[0].Length > 1))
            {
                mincols(ref rowvec, ref res, src, 1.0e+99);
            }
            else if (src[0].Length == 1)
            {
                mincols(ref rowvec, ref res, src, 1.0e+99);
            }
            else
            {
                minrows(ref rowvec, ref res, src, 1.0e+99);
            }
        }

        /// <summary>
        /// 矩阵希尔排序法
        /// </summary>
        /// <param name="arr"></param>
        /// <param name="inx"></param>
        /// <exception cref="Exception"></exception>
        public static void my_shell(ref double[][] arr, ref int[] inx)
        {
            int n = inx.Length;
            if (length_m(arr) != n)
            {
                throw new Exception("MY_SHELL: data and index length don't match");
            }

            var values = new double[n];
            for (int i = 0; i < n; i++)
            {
                values[i] = get(arr, i);
                inx[i] = i;
            }

            Array.Sort(values, inx);

            for (int i = 0; i < n; i++)
            {
                set(ref arr, i, values[i]);
            }
        }

        /// <summary>
        /// 从给定的二维矩阵 matrix 中删除由 delinx 指定的元素，并调整矩阵的大小以反映删除后的结果
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="delinx"></param>
        /// <exception cref="Exception"></exception>
        public static void matselect(ref double[][] matrix, int[] delinx)
        {
            int rows = matrix.Length;
            int cols = matrix[0].Length;
            int total = rows * cols;
            if (delinx.Length > total)
            {
                throw new Exception("DCOLList<double[]>SELECT: more indices to delete than existent !");
            }

            var indices = (int[])delinx.Clone();
            Array.Sort(indices);

            var keep = new List<double>(total - indices.Length);
            int deleteIndex = 0;
            for (int i = 0; i < total; i++)
            {
                if (deleteIndex < indices.Length && indices[deleteIndex] == i)
                {
                    deleteIndex++;
                    continue;
                }

                keep.Add(get(matrix, i));
            }

            if (rows > 1 && cols == 1)
            {
                matrix = new double[0][].Resize(keep.Count, 1);
            }
            else
            {
                matrix = new double[0][].Resize(1, keep.Count);
            }

            for (int i = 0; i < keep.Count; i++)
            {
                set(ref matrix, i, keep[i]);
            }
        }

        /// <summary>
        /// 矩阵重排，去重（有点复杂）
        /// </summary>
        /// <param name="v3"></param>
        /// <param name="v1"></param>
        /// <param name="v2_orig"></param>
        /// <param name="v2"></param>
        /// <param name="ind"></param>
        public static void v2sort_m(ref double[] v3, double[] v1, double[] v2_orig, out double[][] v2, ref int[] ind)
        {
            v2 = new double[0][].Resize(v2_orig.Length, 1);
            copy_m_v(ref v2, v2_orig);

            int lv1 = v1.Length;
            int lv2 = length_m(v2);
            int lv2init = lv2;

            double[][] v1ones = new double[0][].Resize(1, lv1);
            double[][] v2ones = new double[0][].Resize(lv2, 1);
            ones_matrix(ref v1ones);
            ones_matrix(ref v2ones);

            double[][] help3 = new double[0][].Resize(1, lv1);
            copy_m_v(ref help3, v1);
            var help1 = matdotmat(v2ones, help3);
            var help2 = matdotmat(v2, v1ones);
            sub_m_m(ref help1, help2);
            var diff = abs_m_m(help1);

            if (lv1 > lv2)
            {
                double[][] temp = new double[0][].Resize(1, lv1);
                ind = new int[lv1];

                var dum = mul_m_ms(temp, -1.0);
                min_mat(ref temp, ref ind, diff);

                for (int i = 0; i < ind.Length; i++)
                {
                    ind[i] = i;
                }

                my_shell(ref dum, ref ind);

                for (int i = 0; i < lv1 - lv2init; i++)
                {
                    int i2 = Math.Min(lv2, ind[i]);
                    if ((ind[i] <= lv2init) && (i + 1 < lv2init))
                    {
                        double[][] help4 = new double[0][].Resize(i2 + 1, 1);
                        double[][] v1_m = new double[0][].Resize(v1.Length, 1);
                        copy_m_v(ref v1_m, v1);
                        connect_rows(ref help4, v2, 0, i2 - 1, v1_m, ind[i], ind[i]);

                        double[][] help5 = new double[0][].Resize(help4.Length + lv2 - i2, 1);
                        connect_rows(ref help5, help4, 0, help4.Length - 1, v2, i2, lv2 - 1);
                        v2 = help5;
                    }
                    else
                    {
                        double[][] help4 = new double[0][].Resize(v2.Length + 1, 1);
                        double[][] v1_m = new double[0][].Resize(v1.Length, 1);
                        copy_m_v(ref v1_m, v1);
                        connect_rows(ref help4, v2, 0, v2.Length - 1, v1_m, lv2, lv2);
                        v2 = help4;
                    }

                    lv2++;
                    dum = new double[0][].Resize(1, lv2);
                    temp = new double[0][].Resize(1, lv2);
                    ind = new int[lv2];

                    var diffT = transpose(diff);
                    min_mat(ref temp, ref ind, diffT);
                    dum = mul_m_ms(temp, -1.0);

                    for (int j = 0; j < ind.Length; j++)
                    {
                        ind[j] = j;
                    }

                    my_shell(ref dum, ref ind);

                    int trim = Math.Max(0, lv2 - lv1 - 1);
                    var ind_help = new int[trim];
                    Array.Copy(ind, ind_help, trim);
                    Array.Sort(ind_help);
                    matselect(ref v2, ind_help);
                }
            }

            copy_v_m(ref v3, v2);
        }

        /// <summary>
        /// 找到矩阵最近的一个val的索引
        /// </summary>
        /// <param name="v"></param>
        /// <param name="val"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static int find_s(double[][] v, double val)
        {
            if ((v.Length != 1) && (v[0].Length != 1))
            {
                throw new Exception("Input must be a vector !");
            }

            int size = length_m(v);
            for (int i = 0; i < size; i++)
            {
                if (get(v, i) == val)
                {
                    return i;
                }
            }

            return -1;
        }

        /// <summary>
        /// 矩阵列求和
        /// </summary>
        /// <param name="src"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] sumcols(double[][] src)
        {
            var sums = src.ToMatrix().ColumnSums().ToArray();
            var rowvec = new double[0][].Resize(1, sums.Length);
            for (int i = 0; i < sums.Length; i++)
            {
                rowvec[0][i] = sums[i];
            }

            return rowvec;
        }

        /// <summary>
        /// 矩阵行求和
        /// </summary>
        /// <param name="src"></param>
        /// <exception cref="Exception"></exception>
        public static double[][] sumrows(double[][] src)
        {
            var sums = src.ToMatrix().RowSums().ToArray();
            var colvec = new double[0][].Resize(sums.Length, 1);
            for (int i = 0; i < sums.Length; i++)
            {
                colvec[i][0] = sums[i];
            }

            return colvec;
        }

        /// <summary>
        /// 被视为矩阵的数组索引
        /// </summary>
        /// <param name="vec"></param>
        /// <param name="index"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static double get(double[][] vec, int index)
        {
            int cols = vec[0].Length;
            int rows = vec.Length;
            if ((rows > 1) && (cols > 1))
            {
                throw new Exception("must have vector argument!");
            }

            if ((index < 0) || (index >= rows * cols))
            {
                throw new Exception("index exceeds matrix-dimensions!");
            }

            return cols > 1 ? vec[0][index] : vec[index][0];
        }

        /// <summary>
        /// 被视为矩阵的数组赋值
        /// </summary>
        /// <param name="vec"></param>
        /// <param name="index"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static void set(ref double[][] vec, int index, double val)
        {
            int cols = vec[0].Length;
            int rows = vec.Length;
            if ((rows > 1) && (cols > 1))
            {
                throw new Exception("must have vector argument!");
            }

            if ((index < 0) || (index >= rows * cols))
            {
                throw new Exception("index exceeds matrix-dimensions!");
            }

            if (cols > 1)
            {
                vec[0][index] = val;
            }
            else
            {
                vec[index][0] = val;
            }
        }

        /// <summary>
        /// 矩阵向量拷贝，数组被视为矩阵
        /// </summary>
        /// <param name="m"></param>
        /// <param name="v"></param>
        /// <exception cref="Exception"></exception>
        public static void copy_m_v(ref double[][] m, double[] v)
        {
            int cols = m[0].Length;
            int rows = m.Length;
            int size = v.Length;

            if (Math.Min(cols, rows) != 1)
            {
                throw new Exception("matrix must be column or row");
            }

            if (Math.Max(cols, rows) != size)
            {
                throw new Exception("matrix and vector do not match");
            }

            for (int i = 0; i < size; i++)
            {
                int r = cols > 1 ? 0 : i;
                int c = cols > 1 ? i : 0;
                m[r][c] = v[i];
            }
        }

        /// <summary>
        /// 矩阵向量拷贝，数组被视为矩阵
        /// </summary>
        /// <param name="v"></param>
        /// <param name="m"></param>
        /// <exception cref="Exception"></exception>
        public static void copy_v_m(ref double[] v, double[][] m)
        {
            int cols = m[0].Length;
            int rows = m.Length;
            int size = v.Length;
            if (Math.Min(cols, rows) != 1)
            {
                throw new Exception("matrix must be column or row");
            }

            if (Math.Max(cols, rows) != size)
            {
                throw new Exception("matrix and vector do not match");
            }

            for (int i = 0; i < size; i++)
            {
                int r = cols > 1 ? 0 : i;
                int c = cols > 1 ? i : 0;
                v[i] = m[r][c];
            }
        }

        /// <summary>
        /// 获取矩阵对角值
        /// </summary>
        /// <param name="diag"></param>
        /// <param name="mat"></param>
        /// <param name="offset"></param>
        /// <exception cref="Exception"></exception>
        public static void getDiag(ref double[] diag, double[][] mat, int offset = 0)
        {
            var matrix = mat.ToMatrix();
            int size = Math.Min(matrix.RowCount, matrix.ColumnCount) - Math.Abs(offset);
            if (size < 0)
            {
                throw new Exception("offset larger than matrix size");
            }

            diag = matrix.Diagonal(offset).ToArray();
        }

        /// <summary>
        /// 矩阵乘向量
        /// </summary>
        /// <param name="mat"></param>
        /// <param name="input"></param>
        public static double[] matdotvec(double[][] mat, double[] input)
        {
            return mat.ToMatrix().Multiply(input.ToVector()).ToArray();
        }

        /// <summary>
        /// 向量乘矩阵
        /// </summary>
        /// <param name="mat"></param>
        /// <param name="input"></param>
        public static double[] vecdotmat(double[][] mat, double[] input)
        {
            return input.ToVector().ToRowMatrix().Multiply(mat.ToMatrix()).Row(0).ToArray();
        }

        private static double sqrInternal(double value)
        {
            return value * value;
        }

        private static void ReplaceMatrix(ref double[][] target, Matrix<double> matrix)
        {
            target = matrix.ToJaggedArray();
        }
    }
}
