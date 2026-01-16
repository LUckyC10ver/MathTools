using System;
using System.IO;

namespace Math.Core.Polynomials
{
    /// <summary>
    /// 多项式系数与基本运算。
    /// </summary>
    public class Polynomial : PolynomialBase
    {
        private double[] _coeffs;

        /// <summary>
        /// 默认构造函数。
        /// </summary>
        public Polynomial()
        {
            _coeffs = Array.Empty<double>();
        }

        /// <summary>
        /// 以阶数构造多项式。
        /// </summary>
        public Polynomial(int degree)
        {
            if (degree < 0)
            {
                throw new Exception("degree must be non-negative");
            }

            _coeffs = new double[degree];
        }

        /// <inheritdoc />
        public override string ClassName => "Polynomial";

        /// <inheritdoc />
        public override int Degree => _coeffs.Length - 1;

        /// <summary>
        /// 系数数量。
        /// </summary>
        public int Size => _coeffs.Length;

        /// <summary>
        /// 访问系数。
        /// </summary>
        public double this[int index]
        {
            get => _coeffs[index];
            set => _coeffs[index] = value;
        }

        /// <inheritdoc />
        public override void Output(TextWriter writer)
        {
            writer.WriteLine("class Polynomial");
            writer.WriteLine("polynomial coefficients:");
            for (int i = 0; i < _coeffs.Length; i++)
            {
                writer.WriteLine(_coeffs[i]);
            }
        }

        /// <inheritdoc />
        public override double Invoke(double argument)
        {
            double result = 0.0;
            for (int i = _coeffs.Length - 1; i >= 0; i--)
            {
                result = _coeffs[i] + result * argument;
            }

            return result;
        }

        /// <inheritdoc />
        public override UnaryFunction GetDerivativeFunction()
        {
            var derived = Clone();
            derived.Differentiate();
            return derived;
        }

        /// <inheritdoc />
        public override double GetDerivativeValue(double argument)
        {
            var derived = Clone();
            derived.Differentiate();
            return derived.Invoke(argument);
        }

        /// <inheritdoc />
        public override void SetCoefficient(int index, double coeff)
        {
            _coeffs[index] = coeff;
        }

        /// <summary>
        /// 对多项式求导。
        /// </summary>
        public void Differentiate()
        {
            if (_coeffs.Length == 0)
            {
                return;
            }

            int size = _coeffs.Length - 1;
            if (size == 0)
            {
                _coeffs = Array.Empty<double>();
                return;
            }

            var help = new double[size];
            for (int i = 0; i < size; i++)
            {
                help[i] = _coeffs[i + 1] * (i + 1);
            }

            _coeffs = help;
        }

        /// <summary>
        /// 对多项式积分。
        /// </summary>
        public void Integrate()
        {
            int size = _coeffs.Length + 1;
            var help = new double[size];
            help[0] = 0.0;

            for (int i = 1; i < size; i++)
            {
                help[i] = _coeffs[i - 1] / i;
            }

            _coeffs = help;
        }

        /// <summary>
        /// 重置多项式阶数。
        /// </summary>
        public void Resize(int size)
        {
            if (size < 0)
            {
                throw new Exception("size must be non-negative");
            }

            Array.Resize(ref _coeffs, size);
        }

        /// <summary>
        /// 多项式乘法并写回自身。
        /// </summary>
        public void MultiplyAssign(Polynomial rhs)
        {
            int newSize = rhs.Size + Size - 1;
            var newCoeffs = new double[newSize];

            for (int i = 0; i < rhs.Size; i++)
            {
                for (int j = 0; j < Size; j++)
                {
                    int index = i + j;
                    newCoeffs[index] += rhs._coeffs[i] * _coeffs[j];
                }
            }

            _coeffs = newCoeffs;
        }

        /// <summary>
        /// 多项式加法并写回自身。
        /// </summary>
        public void AddAssign(Polynomial rhs)
        {
            int newSize = rhs.Size > Size ? rhs.Size : Size;
            Array.Resize(ref _coeffs, newSize);

            for (int i = 0; i < rhs.Size; i++)
            {
                _coeffs[i] += rhs._coeffs[i];
            }
        }

        /// <summary>
        /// 克隆多项式。
        /// </summary>
        public Polynomial Clone()
        {
            var clone = new Polynomial(0);
            clone._coeffs = new double[_coeffs.Length];
            Array.Copy(_coeffs, clone._coeffs, _coeffs.Length);
            return clone;
        }
    }
}
