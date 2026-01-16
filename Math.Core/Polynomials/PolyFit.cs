using System;
using System.IO;
namespace MathTools.Core.Polynomials
{
    /// <summary>
    /// 多项式最小二乘拟合。
    /// </summary>
    public class PolyFit : DifferentiableFunction
    {
        private readonly double[][] _x;
        private readonly double[] _r;
        private readonly Polynomial _polynom;
        private readonly Polynomial _derivedPolynom;
        private readonly Polynomial _integratedPolynom;
        private readonly bool[] _enabledCoefficients;
        private double _xMin;
        private double _xMax;
        private bool _changed;

        /// <summary>
        /// 构造指定系数数量的拟合器。
        /// </summary>
        public PolyFit(int noCoefficients = 2)
        {
            _x = new double[noCoefficients][];
            for (int i = 0; i < noCoefficients; i++)
            {
                _x[i] = new double[noCoefficients];
            }

            _r = new double[noCoefficients];
            _polynom = new Polynomial(noCoefficients);
            _derivedPolynom = new Polynomial(Math.Max(noCoefficients - 1, 0));
            _integratedPolynom = new Polynomial(noCoefficients + 1);
            _enabledCoefficients = new bool[noCoefficients];
            for (int i = 0; i < _enabledCoefficients.Length; i++)
            {
                _enabledCoefficients[i] = true;
            }
        }

        /// <inheritdoc />
        public override string ClassName => "PolyFit";

        /// <summary>
        /// 添加一个点。
        /// </summary>
        public void AddPoint(double x, double y)
        {
            if (_x[0][0] == 0.0)
            {
                _xMin = x;
                _xMax = x;
            }
            else
            {
                if (x < _xMin)
                {
                    _xMin = x;
                }

                if (x > _xMax)
                {
                    _xMax = x;
                }
            }

            _x[0][0] += 1.0;
            double help = x;

            for (int s = 1; s < _x.Length; s++)
            {
                int row = s;
                int col = 0;
                while (row >= 0)
                {
                    _x[row][col] += help;
                    row--;
                    col++;
                }

                help *= x;
            }

            for (int s = 1; s < _x.Length; s++)
            {
                int row = _x.Length - 1;
                int col = s;
                while (col < _x.Length)
                {
                    _x[row][col] += help;
                    row--;
                    col++;
                }

                help *= x;
            }

            help = y;
            for (int s = 0; s < _x.Length; s++)
            {
                _r[s] += help;
                help *= x;
            }

            SetChanged();
        }

        /// <summary>
        /// 计算拟合多项式。
        /// </summary>
        public void DoFitting()
        {
            if (!_changed)
            {
                return;
            }

            int activeCoefficients = 0;
            var reducedIndex = new int[_enabledCoefficients.Length];
            for (int i = 0; i < reducedIndex.Length; i++)
            {
                reducedIndex[i] = -1;
            }

            for (int i = 0; i < _enabledCoefficients.Length; i++)
            {
                if (_enabledCoefficients[i])
                {
                    reducedIndex[i] = activeCoefficients;
                    activeCoefficients++;
                }
            }

            var lMatrix = new double[activeCoefficients][];
            for (int i = 0; i < activeCoefficients; i++)
            {
                lMatrix[i] = new double[activeCoefficients];
            }

            var indx = new int[activeCoefficients];
            var rhs = new double[activeCoefficients];

            for (int i = 0; i < _enabledCoefficients.Length; i++)
            {
                for (int k = 0; k < _enabledCoefficients.Length; k++)
                {
                    if (reducedIndex[i] != -1 && reducedIndex[k] != -1)
                    {
                        lMatrix[reducedIndex[i]][reducedIndex[k]] = _x[i][k];
                    }
                }
            }

            for (int i = 0; i < _enabledCoefficients.Length; i++)
            {
                if (reducedIndex[i] != -1)
                {
                    rhs[reducedIndex[i]] = _r[i];
                }
            }

            if (_x[0][0] >= _x.Length)
            {
                Functions.LUDecompose(lMatrix, indx, out _);
                Functions.LUBacksubstitute(lMatrix, indx, rhs);

                for (int i = 0; i < _enabledCoefficients.Length; i++)
                {
                    if (reducedIndex[i] != -1)
                    {
                        _polynom[i] = rhs[reducedIndex[i]];
                    }
                    else
                    {
                        _polynom[i] = 0.0;
                    }
                }

                var derived = _polynom.Clone();
                derived.Differentiate();
                CopyPolynomial(_derivedPolynom, derived);

                var integrated = _polynom.Clone();
                integrated.Integrate();
                CopyPolynomial(_integratedPolynom, integrated);

                ResetChanged();
            }
        }

        /// <summary>
        /// 禁用指定系数。
        /// </summary>
        public void DisableCoefficient(int coefficientIndex)
        {
            if (_enabledCoefficients[coefficientIndex])
            {
                _enabledCoefficients[coefficientIndex] = false;
                SetChanged();
            }
        }

        /// <summary>
        /// 启用指定系数。
        /// </summary>
        public void EnableCoefficient(int coefficientIndex)
        {
            if (!_enabledCoefficients[coefficientIndex])
            {
                _enabledCoefficients[coefficientIndex] = true;
                SetChanged();
            }
        }

        /// <summary>
        /// 复位内部统计量。
        /// </summary>
        public void Reset()
        {
            for (int i = 0; i < _x.Length; i++)
            {
                Array.Clear(_x[i], 0, _x[i].Length);
                _r[i] = 0.0;
            }

            _xMin = 0.0;
            _xMax = 0.0;
        }

        /// <summary>
        /// 获取最小 x 值。
        /// </summary>
        public double XMin => _xMin;

        /// <summary>
        /// 获取最大 x 值。
        /// </summary>
        public double XMax => _xMax;

        /// <summary>
        /// 获取点数量。
        /// </summary>
        public int GetNoPoints()
        {
            return (int)_x[0][0];
        }

        /// <summary>
        /// 获取拟合多项式。
        /// </summary>
        public Polynomial GetPolynomial()
        {
            DoFitting();
            return _polynom.Clone();
        }

        /// <summary>
        /// 获取拟合多项式引用。
        /// </summary>
        public Polynomial GetPolynomialRef()
        {
            return _polynom;
        }

        /// <summary>
        /// 获取导数多项式。
        /// </summary>
        public Polynomial GetIntegratedPolynomial()
        {
            DoFitting();
            return _integratedPolynom.Clone();
        }

        /// <inheritdoc />
        public override UnaryFunction GetDerivativeFunction()
        {
            return _derivedPolynom;
        }

        /// <inheritdoc />
        public override double Invoke(double argument)
        {
            return _polynom.Invoke(argument);
        }

        /// <inheritdoc />
        public override void Output(TextWriter writer)
        {
            writer.WriteLine("X =");
            for (int i = 0; i < _x.Length; i++)
            {
                for (int j = 0; j < _x[i].Length; j++)
                {
                    writer.Write($"{_x[i][j]} ");
                }

                writer.WriteLine();
            }

            writer.WriteLine("r =");
            for (int i = 0; i < _r.Length; i++)
            {
                writer.WriteLine(_r[i]);
            }

            writer.WriteLine("polynom =");
            _polynom.Output(writer);
        }

        /// <summary>
        /// 获取导数值。
        /// </summary>
        public override double GetDerivative(double argument)
        {
            return _derivedPolynom.Invoke(argument);
        }

        /// <summary>
        /// 计算积分。
        /// </summary>
        public double Integral(double a, double b)
        {
            DoFitting();
            return _integratedPolynom.Invoke(b) - _integratedPolynom.Invoke(a);
        }

        private void ResetChanged()
        {
            _changed = false;
        }

        private void SetChanged()
        {
            _changed = true;
        }

        private static void CopyPolynomial(Polynomial target, Polynomial source)
        {
            target.Resize(source.Size);
            for (int i = 0; i < source.Size; i++)
            {
                target[i] = source[i];
            }
        }
    }
}
