using System;

namespace Math.Core.Legacy.Regression
{
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

        /// <summary>
        /// 添加样本点 (x,y,z)。
        /// </summary>
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

        /// <summary>
        /// 更新回归系数。
        /// </summary>
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

        /// <summary>
        /// 获取回归系数 (a0, a1, a2)。
        /// </summary>
        public double[] Coefficients => (double[])_coefficients.Clone();

        /// <summary>
        /// 设定系数启用状态。
        /// </summary>
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

        /// <summary>
        /// 重置状态。
        /// </summary>
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
                _coefficients[0] = 0.0;
                _coefficients[1] = 0.0;
                _coefficients[2] = 0.0;
                _coefficients[indices[0]] = sub[0];
            }
        }
    }
}
