using System;

namespace MathTools.Core
{
    public partial class Functions
    {
        /// <summary>
        /// 乘方运算
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static double pow_new(double x, double y)
        {
            if ((x == 0) && (y == 0))
            {
                return 1;
            }

            return Math.Pow(x, y);
        }

        /// <summary>
        /// 平方运算
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double sqr(double x)
        {
            return x * x;
        }

        /// <summary>
        /// 交换两个双浮点数
        /// </summary>
        /// <param name="x1"></param>
        /// <param name="x2"></param>
        public static void flip<T>(ref T x1, ref T x2)
        {
            (x2, x1) = (x1, x2);
        }
    }
}
