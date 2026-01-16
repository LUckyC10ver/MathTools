using System;
using System.IO;

namespace Math.Core.Legacy.Polynomials
{
    /// <summary>
    /// 固定最大阶数的多项式（按降阶存储系数）。
    /// </summary>
    public sealed class BCTinyPolynomial
    {
        private readonly double[] _coeff;
        private int _degree;

        /// <summary>
        /// 创建一个固定阶数多项式，系数默认全为零。
        /// </summary>
        public BCTinyPolynomial(int maxDegree)
        {
            if (maxDegree < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(maxDegree));
            }

            _coeff = new double[maxDegree + 1];
            ComputeDegree();
        }

        /// <summary>
        /// 使用给定系数初始化，按降阶存储。
        /// </summary>
        public BCTinyPolynomial(double[] coefficients)
        {
            if (coefficients == null || coefficients.Length == 0)
            {
                throw new ArgumentException("coefficients required", nameof(coefficients));
            }

            _coeff = (double[])coefficients.Clone();
            ComputeDegree();
        }

        /// <summary>
        /// 最大阶数。
        /// </summary>
        public int MaxDegree => _coeff.Length - 1;

        /// <summary>
        /// 当前阶数（忽略最高次的零系数）。
        /// </summary>
        public int Degree => _degree;

        /// <summary>
        /// 获取指定幂次的系数。
        /// </summary>
        public double GetCoefficient(int degree)
        {
            if (degree < 0 || degree > MaxDegree)
            {
                throw new ArgumentOutOfRangeException(nameof(degree));
            }

            return _coeff[MaxDegree - degree];
        }

        /// <summary>
        /// 设置指定幂次的系数。
        /// </summary>
        public void SetCoefficient(int degree, double value)
        {
            if (degree < 0 || degree > MaxDegree)
            {
                throw new ArgumentOutOfRangeException(nameof(degree));
            }

            _coeff[MaxDegree - degree] = value;
            ComputeDegree();
        }

        /// <summary>
        /// 计算多项式值。
        /// </summary>
        public double Evaluate(double argument)
        {
            double result = 0.0;
            for (int i = 0; i < _coeff.Length; i++)
            {
                result = result * argument + _coeff[i];
            }

            return result;
        }

        /// <summary>
        /// 计算导数值。
        /// </summary>
        public double EvaluateDerivative(double argument)
        {
            double result = 0.0;
            for (int i = 0; i < MaxDegree; i++)
            {
                int power = MaxDegree - i;
                result = result * argument + power * _coeff[i];
            }

            return result;
        }

        /// <summary>
        /// 获取导数多项式。
        /// </summary>
        public BCTinyPolynomial GetDerivative()
        {
            if (MaxDegree == 0)
            {
                return new BCTinyPolynomial(new[] { 0.0 });
            }

            var coeff = new double[MaxDegree];
            for (int i = 0; i < MaxDegree; i++)
            {
                int power = MaxDegree - i;
                coeff[i] = power * _coeff[i];
            }

            return new BCTinyPolynomial(coeff);
        }

        /// <summary>
        /// 输出多项式信息。
        /// </summary>
        public void Output(TextWriter writer)
        {
            if (writer == null)
            {
                throw new ArgumentNullException(nameof(writer));
            }

            writer.WriteLine("BCTinyPolynomial");
            writer.WriteLine($"Degree: {Degree}");
        }

        private void ComputeDegree()
        {
            _degree = -1;
            for (int i = 0; i < _coeff.Length; i++)
            {
                if (Math.Abs(_coeff[i]) > 0.0)
                {
                    _degree = MaxDegree - i;
                    return;
                }
            }
        }
    }
}
