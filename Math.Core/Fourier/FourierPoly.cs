using System;
using System.Numerics;

namespace MathTools.Core.Fourier
{
    /// <summary>
    /// 多项式脉冲的傅里叶系数。
    /// </summary>
    public class FourierPoly : FourierSeries
    {
        /// <summary>
        /// 默认构造函数。
        /// </summary>
        public FourierPoly()
        {
        }

        /// <summary>
        /// 计算多项式脉冲的傅里叶系数。
        /// </summary>
        /// <param name="ino">系数数量。</param>
        /// <param name="iorder">多项式阶数（Nullordnung）。</param>
        /// <param name="dmax">半脉冲宽度。</param>
        /// <param name="dshift">中心位移。</param>
        public FourierPoly(int ino, int iorder, double dmax, double dshift)
            : base(ino, Complex.Zero)
        {
            if (dmax == 0.0 || ino == 0 || iorder <= 0)
            {
                return;
            }

            double[,] m = new double[iorder, iorder];
            for (int k1 = 0; k1 < iorder; k1++)
            {
                m[0, k1] = 1.0;
                m[k1, 0] = 1.0;
            }

            for (int row1 = 1; row1 < iorder; row1++)
            {
                for (int col = 1; col < iorder - row1; col++)
                {
                    m[row1, col] = m[row1, col - 1] + m[row1 - 1, col];
                }
            }

            int porder = iorder * 2;
            var poly = new Polynomial(porder);
            double sign = 1.0;
            for (int row2 = 0; row2 < iorder; row2++)
            {
                poly[2 * row2] = sign * m[row2, iorder - 1 - row2];
                sign *= -1.0;
            }

            var poly0 = poly.Clone();
            poly0.Integrate();

            double rectNorm = 1.0 / poly0.Evaluate(1.0);
            m_data[0] = dmax / MathCoreInfo.Pi;

            var fa = new double[porder];
            double koeff = 1.0;
            double norm = 1.0 / dmax;
            for (int k2 = 0; k2 < porder; k2++)
            {
                fa[k2] = poly.Evaluate(1.0) * koeff;
                poly.Differentiate();
                koeff *= norm;
            }

            for (int n = 1; n < Size; n++)
            {
                double dn = n;
                double multiplikator = -1.0 / (dn * dn);
                double koeffizient = 1.0;
                double sg = 0.0;
                double su = 0.0;

                for (int k2 = 0; k2 < iorder; k2++)
                {
                    sg += koeffizient * fa[2 * k2];
                    su += koeffizient * fa[2 * k2 + 1];
                    koeffizient *= multiplikator;
                }

                m_data[n] =
                    -(dn * Math.Sin(dn * dmax) * sg + Math.Cos(dn * dmax) * su)
                    * rectNorm * multiplikator / MathCoreInfo.Pi;
            }

            Move(dshift);
        }

        private sealed class Polynomial
        {
            private readonly double[] _coeffs;

            public Polynomial(int order)
            {
                _coeffs = new double[order];
            }

            public double this[int index]
            {
                get => _coeffs[index];
                set => _coeffs[index] = value;
            }

            public Polynomial Clone()
            {
                var clone = new Polynomial(_coeffs.Length);
                Array.Copy(_coeffs, clone._coeffs, _coeffs.Length);
                return clone;
            }

            public double Evaluate(double x)
            {
                double sum = 0.0;
                double pow = 1.0;
                for (int i = 0; i < _coeffs.Length; i++)
                {
                    sum += _coeffs[i] * pow;
                    pow *= x;
                }

                return sum;
            }

            public void Differentiate()
            {
                for (int i = 0; i < _coeffs.Length - 1; i++)
                {
                    _coeffs[i] = _coeffs[i + 1] * (i + 1);
                }

                _coeffs[^1] = 0.0;
            }

            public void Integrate()
            {
                for (int i = _coeffs.Length - 1; i > 0; i--)
                {
                    _coeffs[i] = _coeffs[i - 1] / i;
                }

                _coeffs[0] = 0.0;
            }
        }
    }
}
