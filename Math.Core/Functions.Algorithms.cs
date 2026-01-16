using System;
using System.Linq;
using Math.Core.Extensions;

namespace Math.Core
{
    public partial class Functions
    {
        /// <summary>
        /// 二分查找，返回找到的索引，不存在返回 -1。
        /// </summary>
        public static int binary_find(double[] values, double value)
        {
            if (values == null || values.Length == 0)
            {
                return -1;
            }

            int index = Array.BinarySearch(values, value);
            return index >= 0 ? index : -1;
        }

        /// <summary>
        /// 向量距离（欧氏距离）。
        /// </summary>
        public static double dist(double[] x, double[] y)
        {
            if (x == null || y == null || x.Length != y.Length)
            {
                throw new Exception("vector sizes differ in dist(BCVector&,BCVector&)");
            }

            return x.ToVector().Subtract(y.ToVector()).L2Norm();
        }

        /// <summary>
        /// 求向量元素和。
        /// </summary>
        public static double suma(double[] a)
        {
            if (a == null || a.Length == 0)
            {
                return 0.0;
            }

            return a.Sum();
        }

        /// <summary>
        /// 获取向量最小值。
        /// </summary>
        public static double min(double[] vec)
        {
            if (vec == null || vec.Length == 0)
            {
                return 0.0;
            }

            double result = vec[0];
            for (int i = 1; i < vec.Length; i++)
            {
                if (vec[i] < result)
                {
                    result = vec[i];
                }
            }

            return result;
        }

        /// <summary>
        /// 获取向量最大值。
        /// </summary>
        public static double max(double[] vec)
        {
            if (vec == null || vec.Length == 0)
            {
                return 0.0;
            }

            double result = vec[0];
            for (int i = 1; i < vec.Length; i++)
            {
                if (vec[i] > result)
                {
                    result = vec[i];
                }
            }

            return result;
        }

        /// <summary>
        /// 插入排序（适合小规模数组）。
        /// </summary>
        public static void piksrt(double[] values)
        {
            if (values == null || values.Length < 2)
            {
                return;
            }

            for (int j = 1; j < values.Length; j++)
            {
                double a = values[j];
                int i = j - 1;
                while (i >= 0 && a < values[i])
                {
                    values[i + 1] = values[i];
                    i--;
                }

                values[i + 1] = a;
            }
        }
    }
}
