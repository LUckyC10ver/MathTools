using System;

namespace Math.Core.Legacy.Regression
{
    /// <summary>
    /// 1~3 维线性方程求解器。
    /// </summary>
    public static class BCTinyLinEquations
    {
        public static double[] Solve1(double[,] coeff, double[] rhs)
        {
            if (coeff[0, 0] == 0.0)
            {
                throw new InvalidOperationException("coeffMatrix(0,0) is zero");
            }

            return new[] { rhs[0] / coeff[0, 0] };
        }

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

    /// <summary>
    /// 一元线性回归：拟合 y = a0 + a1*x。
    /// </summary>
    public sealed class BCTinyLinReg
    {
        private double _sumX;
        private double _sumX2;
        private double _sumY;
        private double _sumYX;
        private int _count;
        private double _a0;
        private double _a1;
        private bool _changed;

        public void AddPoint(double x, double y)
        {
            _sumX += x;
            _sumX2 += x * x;
            _sumY += y;
            _sumYX += y * x;
            _count++;
            _changed = true;
        }

        public void Calculate()
        {
            if (!_changed || _count == 0)
            {
                return;
            }

            double denom = _count * _sumX2 - _sumX * _sumX;
            if (Math.Abs(denom) > 1e-20)
            {
                _a1 = (_count * _sumYX - _sumY * _sumX) / denom;
                _a0 = (_sumY - _a1 * _sumX) / _count;
                _changed = false;
            }
        }

        public double A0 => _a0;
        public double A1 => _a1;

        public void Reset()
        {
            _sumX = 0.0;
            _sumX2 = 0.0;
            _sumY = 0.0;
            _sumYX = 0.0;
            _count = 0;
            _a0 = 0.0;
            _a1 = 0.0;
            _changed = false;
        }
    }

    /// <summary>
    /// 二元线性回归：拟合 z = a0 + a1*x + a2*y。
    /// </summary>
    public sealed class BCTinyLinReg2
    {
        private readonly double[,] _coeffMatrix = new double[3, 3];
        private readonly double[] _rhs = new double[3];
        private readonly bool[] _coefficientStatus = { true, true, true };
        private readonly double[] _coefficients = new double[3];
        private bool _changed;

        public void AddPoint(double x, double y, double z)
        {
            _coeffMatrix[0, 0] += 1.0;
            _rhs[0] += z;

            _coeffMatrix[0, 1] += x;
            _coeffMatrix[1, 0] += x;
            _rhs[1] += z * x;

            _coeffMatrix[0, 2] += y;
            _coeffMatrix[2, 0] += y;
            _rhs[2] += z * y;

            _coeffMatrix[1, 1] += x * x;
            _coeffMatrix[1, 2] += x * y;
            _coeffMatrix[2, 1] += x * y;
            _coeffMatrix[2, 2] += y * y;

            _changed = true;
        }

        public void Calculate()
        {
            if (!_changed || _coeffMatrix[0, 0] == 0.0)
            {
                return;
            }

            int active = 0;
            for (int i = 0; i < _coefficientStatus.Length; i++)
            {
                if (_coefficientStatus[i])
                {
                    active++;
                }
            }

            if (active == 3)
            {
                var solved = BCTinyLinEquations.Solve3(_coeffMatrix, _rhs);
                Array.Copy(solved, _coefficients, 3);
            }
            else if (active == 2)
            {
                SolveReduced(2);
            }
            else if (active == 1)
            {
                SolveReduced(1);
            }
            else
            {
                throw new InvalidOperationException("No active coefficients");
            }

            _changed = false;
        }

        public double[] Coefficients => (double[])_coefficients.Clone();

        public void SetCoefficientStatus(bool a0, bool a1, bool a2)
        {
            if (_coefficientStatus[0] != a0 || _coefficientStatus[1] != a1 || _coefficientStatus[2] != a2)
            {
                _coefficientStatus[0] = a0;
                _coefficientStatus[1] = a1;
                _coefficientStatus[2] = a2;
                _changed = true;
            }
        }

        public void Reset()
        {
            Array.Clear(_coeffMatrix, 0, _coeffMatrix.Length);
            Array.Clear(_rhs, 0, _rhs.Length);
            Array.Clear(_coefficients, 0, _coefficients.Length);
            _changed = false;
        }

        private void SolveReduced(int activeCount)
        {
            int[] indices = new int[activeCount];
            int idx = 0;
            for (int i = 0; i < _coefficientStatus.Length; i++)
            {
                if (_coefficientStatus[i])
                {
                    indices[idx++] = i;
                }
            }

            if (activeCount == 2)
            {
                double[,] subCoeff =
                {
                    { _coeffMatrix[indices[0], indices[0]], _coeffMatrix[indices[0], indices[1]] },
                    { _coeffMatrix[indices[1], indices[0]], _coeffMatrix[indices[1], indices[1]] }
                };
                double[] subRhs = { _rhs[indices[0]], _rhs[indices[1]] };
                double[] sub = BCTinyLinEquations.Solve2(subCoeff, subRhs);
                _coefficients[indices[0]] = sub[0];
                _coefficients[indices[1]] = sub[1];
            }
            else
            {
                double[,] subCoeff = { { _coeffMatrix[indices[0], indices[0]] } };
                double[] subRhs = { _rhs[indices[0]] };
                double[] sub = BCTinyLinEquations.Solve1(subCoeff, subRhs);
                Array.Clear(_coefficients, 0, _coefficients.Length);
                _coefficients[indices[0]] = sub[0];
            }
        }
    }
}
