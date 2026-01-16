using System;
using System.Numerics;

namespace MathTools.Core.Fourier
{
    /// <summary>
    /// 矩形脉冲的傅里叶系数。
    /// </summary>
    public class FourierRect : FourierSeries
    {
        /// <summary>
        /// 构造矩形脉冲的傅里叶系数。
        /// </summary>
        /// <param name="ino">系数数量。</param>
        /// <param name="thetaMin">最小角度（0~2π）。</param>
        /// <param name="thetaMax">最大角度（0~2π）。</param>
        /// <param name="h">脉冲高度。</param>
        public FourierRect(int ino = 0, double thetaMin = 0.0, double thetaMax = 0.0, double h = 1.0)
            : base(ino, Complex.Zero)
        {
            if (ino == 0)
            {
                return;
            }

            EnsureRange(thetaMin, 0.0, 2.0 * MathCoreInfo.Pi, "minimum angle");
            EnsureRange(thetaMax, 0.0, 2.0 * MathCoreInfo.Pi, "maximum angle");

            Complex i = new Complex(0.0, 1.0);

            m_data[0] = (thetaMax - thetaMin) / (2.0 * MathCoreInfo.Pi);
            if (thetaMax < thetaMin)
            {
                m_data[0] += 1.0;
            }

            m_data[0] *= h;

            for (int n = 1; n < Size; n++)
            {
                m_data[n] = h * i / (2.0 * MathCoreInfo.Pi * n)
                            * (Complex.Exp(-i * n * thetaMax) - Complex.Exp(-i * n * thetaMin));
            }
        }

        private static void EnsureRange(double value, double min, double max, string name)
        {
            if (value < min || value > max)
            {
                throw new Exception($"{name} out of range");
            }
        }
    }
}
