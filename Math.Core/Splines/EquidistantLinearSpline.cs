using System;

namespace MathTools.Core.Splines
{
    /// <summary>
    /// 等距线性样条。
    /// </summary>
    public class EquidistantLinearSpline
    {
        private double[] _data = Array.Empty<double>();
        private double _xmin;
        private double _xmax;
        private double _dx = 1.0;
        private double _epsilon = MathCoreInfo.BcEpsEquidistantLinearSpline;

        /// <summary>
        /// 默认构造函数。
        /// </summary>
        public EquidistantLinearSpline()
        {
        }

        /// <summary>
        /// 通过数据区间初始化。
        /// </summary>
        public EquidistantLinearSpline(double[] data, double xmin, double xmax)
        {
            Init(data, 0, data?.Length ?? 0, xmin, xmax);
        }

        /// <summary>
        /// 通过点和值初始化。
        /// </summary>
        public EquidistantLinearSpline(double[] xValues, double[] yValues, bool sorted = false)
        {
            if (xValues == null || yValues == null || xValues.Length != yValues.Length)
            {
                throw new Exception("x values not equidistant");
            }

            var x = (double[])xValues.Clone();
            var y = (double[])yValues.Clone();
            if (!sorted)
            {
                Functions.sortPoints(ref x, ref y);
            }

            if (!CheckEquidistance(x))
            {
                throw new Exception("x values not equidistant");
            }

            _xmin = x[0];
            _xmax = x[^1];
            _dx = (_xmax - _xmin) / (x.Length - 1);
            _epsilon = MathCoreInfo.BcEpsEquidistantLinearSpline * Math.Max(Math.Abs(_xmin), Math.Abs(_xmax));
            _data = y;
        }

        /// <summary>
        /// 初始化数据区间（使用 [start, start+count)）。
        /// </summary>
        public void Init(double[] data, int start, int count, double xmin, double xmax)
        {
            if (data == null)
            {
                throw new Exception("data is null");
            }

            if (xmin > xmax)
            {
                throw new Exception($"invalid range [xmin,xmax]=[{xmin},{xmax}] specified");
            }

            if (count < 1 || start < 0 || start + count > data.Length)
            {
                throw new Exception($"invalid iterator range specified (size={count})");
            }

            _xmin = xmin;
            _xmax = xmax;
            _dx = (_xmax - _xmin) / (count - 1);
            _epsilon = MathCoreInfo.BcEpsEquidistantLinearSpline * Math.Max(Math.Abs(_xmin), Math.Abs(_xmax));

            _data = new double[count];
            Array.Copy(data, start, _data, 0, count);
        }

        /// <summary>
        /// 计算样条值。
        /// </summary>
        public double Evaluate(double x)
        {
            if (x < _xmin - _epsilon || x > _xmax + _epsilon)
            {
                throw new Exception($"input value x={x} out of range [{_xmin},{_xmax}]");
            }

            int index = (int)((x - _xmin) / _dx);
            if (index <= -1)
            {
                return _data[0];
            }

            if (_data.Length - 1 <= index)
            {
                return _data[^1];
            }

            double xLeft = _xmin + index * _dx;
            double dyDx = (_data[index + 1] - _data[index]) / _dx;
            return _data[index] + dyDx * (x - xLeft);
        }

        /// <summary>
        /// 计算积分。
        /// </summary>
        public double GetIntegral(double a, double b)
        {
            if (a < _xmin || a > _xmax)
            {
                throw new Exception($"input value a={a} out of range [{_xmin},{_xmax}]");
            }

            if (b < _xmin || b > _xmax)
            {
                throw new Exception($"input value b={b} out of range [{_xmin},{_xmax}]");
            }

            if (_data.Length < 2)
            {
                throw new Exception($"value size less 2 ({_data.Length})");
            }

            double signum = 1.0;
            if (b < a)
            {
                (a, b) = (b, a);
                signum = -1.0;
            }

            double minIndex = (a - _xmin - _dx / 2.0) / _dx;
            double maxIndex = (b - _xmin - _dx / 2.0) / _dx;
            int indexOfa = minIndex <= 0.0 ? 0 : (int)Math.Floor(minIndex);
            int indexOfb = maxIndex <= 0.0 ? 0 : (int)Math.Floor(maxIndex);

            double result;
            if (indexOfa == indexOfb)
            {
                double xLeft = _xmin + indexOfa * _dx;
                double dyDx = (_data[indexOfa + 1] - _data[indexOfa]) / _dx;
                double ya = _data[indexOfa] + dyDx * (a - xLeft);
                double yb = _data[indexOfa] + dyDx * (b - xLeft);
                result = 0.5 * (ya + yb) * (b - a);
            }
            else
            {
                result = 0.0;
                double xa = _xmin + indexOfa * _dx;
                double xb = _xmin + indexOfb * _dx;

                double dyDx = (_data[indexOfa + 1] - _data[indexOfa]) / _dx;
                double y0 = _data[indexOfa] + dyDx * (a - xa);
                result += 0.5 * (y0 + _data[indexOfa + 1]) * (xa + _dx - a);

                if (indexOfb < _data.Length - 1)
                {
                    dyDx = (_data[indexOfb + 1] - _data[indexOfb]) / _dx;
                }
                else
                {
                    dyDx = (_data[indexOfb] - _data[indexOfb - 1]) / _dx;
                }

                y0 = _data[indexOfb] + dyDx * (b - xb);
                result += 0.5 * (y0 + _data[indexOfb]) * (b - xb);

                if ((indexOfb - indexOfa) > 1)
                {
                    double sum = 0.5 * (_data[indexOfa + 1] + _data[indexOfb]);
                    for (int i = indexOfa + 2; i < indexOfb; i++)
                    {
                        sum += _data[i];
                    }

                    result += sum * _dx;
                }
            }

            return result * signum;
        }

        private bool CheckEquidistance(double[] vec)
        {
            if (vec.Length < 2)
            {
                return false;
            }

            if (vec.Length == 2)
            {
                return true;
            }

            double delta = vec[1] - vec[0];
            for (int i = 2; i < vec.Length; i++)
            {
                double help = vec[i] - vec[i - 1];
                if (Math.Abs(help - delta) > 1e-10)
                {
                    return false;
                }
            }

            return true;
        }
    }
}
