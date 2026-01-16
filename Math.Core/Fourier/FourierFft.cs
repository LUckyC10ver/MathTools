using System;
using System.Numerics;

namespace Math.Core.Fourier
{
    /// <summary>
    /// 使用离散数据序列计算傅里叶系数。
    /// </summary>
    public class FourierFft : FourierSeries
    {
        /// <summary>
        /// 默认构造函数。
        /// </summary>
        public FourierFft()
        {
        }

        /// <summary>
        /// 使用数据序列初始化系数。
        /// </summary>
        /// <param name="ino">系数数量。</param>
        /// <param name="data">数据序列。</param>
        public FourierFft(int ino, double[] data)
            : base(ino, Complex.Zero)
        {
            if (ino == 0)
            {
                return;
            }

            if (data == null || data.Length == 0)
            {
                throw new Exception("Could not initialize without data");
            }

            double deltaTheta = 2.0 * MathCoreInfo.Pi / data.Length;

            for (int angle1 = 0; angle1 < data.Length; angle1++)
            {
                m_data[0] += data[angle1];
            }

            m_data[0] /= data.Length;

            Complex e1 = Complex.Exp(new Complex(0.0, -1.0) * deltaTheta);
            Complex en = Complex.One;

            for (int n = 1; n < ino; n++)
            {
                Complex ek = Complex.One;
                Complex ekPrev = Complex.One;

                en *= e1;

                for (int angle2 = 0; angle2 < data.Length; angle2++)
                {
                    ekPrev = ek;
                    ek *= en;
                    if (data[angle2] != 0.0)
                    {
                        m_data[n] += data[angle2] * (ek - ekPrev);
                    }
                }

                m_data[n] *= new Complex(0.0, 1.0) / (2.0 * MathCoreInfo.Pi * n);
            }
        }
    }
}
