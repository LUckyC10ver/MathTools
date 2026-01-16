using System;

namespace Math.Core.Legacy.Regression
{
    /// <summary>
    /// 1~3 维线性方程求解器。
    /// </summary>
    public static class BCTinyLinEquations
    {
        /// <summary>
        /// 解一元线性方程。
        /// </summary>
        public static double[] Solve1(double[,] coeff, double[] rhs)
        {
            if (coeff[0, 0] == 0.0)
            {
                throw new InvalidOperationException("coeffMatrix(0,0) is zero");
            }

            return new[] { rhs[0] / coeff[0, 0] };
        }

        /// <summary>
        /// 解二元线性方程。
        /// </summary>
        public static double[] Solve2(double[,] coeff, double[] rhs)
        {
            if (coeff[0, 0] == 0.0)
            {
                throw new InvalidOperationException("coeffMatrix(0,0) is zero");
            }

            double subCoeff = coeff[1, 1] * coeff[0, 0] - coeff[1, 0] * coeff[0, 1];
            double subRhs = rhs[1] * coeff[0, 0] - rhs[0] * coeff[1, 0];
            double x1 = subRhs / subCoeff;
            double x0 = (rhs[0] - x1 * coeff[0, 1]) / coeff[0, 0];
            return new[] { x0, x1 };
        }

        /// <summary>
        /// 解三元线性方程。
        /// </summary>
        public static double[] Solve3(double[,] coeff, double[] rhs)
        {
            if (coeff[0, 0] == 0.0)
            {
                throw new InvalidOperationException("coeffMatrix(0,0) is zero");
            }

            double[,] subCoeff =
            {
                { coeff[1, 1] * coeff[0, 0] - coeff[1, 0] * coeff[0, 1], coeff[1, 2] * coeff[0, 0] - coeff[0, 2] * coeff[1, 0] },
                { coeff[2, 1] * coeff[0, 0] - coeff[0, 1] * coeff[2, 0], coeff[2, 2] * coeff[0, 0] - coeff[0, 2] * coeff[2, 0] }
            };
            double[] subRhs =
            {
                rhs[1] * coeff[0, 0] - rhs[0] * coeff[1, 0],
                rhs[2] * coeff[0, 0] - rhs[0] * coeff[2, 0]
            };

            double[] subResult = Solve2(subCoeff, subRhs);
            double x2 = subResult[1];
            double x1 = subResult[0];
            double x0 = (rhs[0] - x1 * coeff[0, 1] - x2 * coeff[0, 2]) / coeff[0, 0];
            return new[] { x0, x1, x2 };
        }
    }
}
