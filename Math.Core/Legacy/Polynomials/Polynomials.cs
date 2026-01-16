using System;
using System.IO;

namespace MathTools.Core
{
    /// <summary>
    /// 固定最大阶数的多项式（按降阶存储系数）。
    /// </summary>
    public sealed class BCTinyPolynomial
    {
        private readonly double[] _coeff;
        private int _degree;

        public BCTinyPolynomial(int maxDegree)
        {
            if (maxDegree < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(maxDegree));
            }

            _coeff = new double[maxDegree + 1];
            ComputeDegree();
        }

        public BCTinyPolynomial(double[] coefficients)
        {
            if (coefficients == null || coefficients.Length == 0)
            {
                throw new ArgumentException("coefficients required", nameof(coefficients));
            }

            _coeff = (double[])coefficients.Clone();
            ComputeDegree();
        }

        public int MaxDegree => _coeff.Length - 1;
        public int Degree => _degree;

        public double GetCoefficient(int degree)
        {
            if (degree < 0 || degree > MaxDegree)
            {
                throw new ArgumentOutOfRangeException(nameof(degree));
            }

            return _coeff[MaxDegree - degree];
        }

        public void SetCoefficient(int degree, double value)
        {
            if (degree < 0 || degree > MaxDegree)
            {
                throw new ArgumentOutOfRangeException(nameof(degree));
            }

            _coeff[MaxDegree - degree] = value;
            ComputeDegree();
        }

        public double Evaluate(double argument)
        {
            double result = 0.0;
            for (int i = 0; i < _coeff.Length; i++)
            {
                result = result * argument + _coeff[i];
            }

            return result;
        }

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
