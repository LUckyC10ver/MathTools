using System;

namespace MathTools.Core
{
    /// <summary>
    /// Gauss-Lobatto 节点线性样条。
    /// </summary>
    public class GaussLobattoSpline
    {
        private double[] _data = Array.Empty<double>();
        private double _xmin;
        private double _xmax;
        private int _mode;

        /// <summary>
        /// 默认构造函数。
        /// </summary>
        public GaussLobattoSpline(int mode = -1)
        {
            _mode = mode > -1 ? mode : 0;
        }

        /// <summary>
        /// 使用区间数据初始化（使用 [start, start+count)）。
        /// </summary>
        public GaussLobattoSpline(double[] data, int start, int count, double xmin, double xmax, int mode = -1)
        {
            Init(data, start, count, xmin, xmax, mode);
        }

        /// <summary>
        /// 使用点和值初始化。
        /// </summary>
        public GaussLobattoSpline(double[] xValues, double[] yValues, bool sorted = false, int mode = -1)
        {
            if (xValues == null || yValues == null || xValues.Length != yValues.Length)
            {
                throw new Exception("vector sizes do not fit");
            }

            _mode = mode > -1 ? mode : 0;

            var x = (double[])xValues.Clone();
            var y = (double[])yValues.Clone();
            if (!sorted)
            {
                Functions.sortPoints(ref x, ref y);
            }

            _xmin = x[0];
            _xmax = x[x.Length - 1];

            _data = new double[x.Length];
            int externIndex = 0;
            double xLeft = x[externIndex];
            double xRight = x[externIndex + 1];
            double dyDx = (y[externIndex + 1] - y[externIndex]) / (xRight - xLeft);

            for (int i = 0; i < _data.Length; i++)
            {
                double xi = XOfIndex(i);
                bool intervalChanged = false;
                for (; xi > x[externIndex + 1] && externIndex < _data.Length - 2; externIndex++)
                {
                    intervalChanged = true;
                }

                if (intervalChanged)
                {
                    xLeft = x[externIndex];
                    xRight = x[externIndex + 1];
                    dyDx = (y[externIndex + 1] - y[externIndex]) / (xRight - xLeft);
                }

                _data[i] = y[externIndex] + dyDx * (xi - x[externIndex]);
            }
        }

        /// <summary>
        /// 返回数据长度。
        /// </summary>
        public int Size => _data.Length;

        /// <summary>
        /// 返回区间最小值。
        /// </summary>
        public double XMin() => _xmin;

        /// <summary>
        /// 返回区间最大值。
        /// </summary>
        public double XMax() => _xmax;

        /// <summary>
        /// 模式参数。
        /// </summary>
        public int Mode => _mode;

        /// <summary>
        /// 初始化数据区间（使用 [start, start+count)）。
        /// </summary>
        public void Init(double[] data, int start, int count, double xmin, double xmax, int mode = -1)
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

            if (mode > -1)
            {
                _mode = mode;
            }

            _xmin = xmin;
            _xmax = xmax;
            _data = new double[count];
            Array.Copy(data, start, _data, 0, count);
        }

        /// <summary>
        /// 计算样条值。
        /// </summary>
        public double Evaluate(double x)
        {
            int index = (int)IndexOfX(x);
            if (index <= -1)
            {
                return _data[0];
            }

            if (index > _data.Length - 1)
            {
                return _data[_data.Length - 1];
            }

            double xLeft = XOfIndex(index);
            double xRight = XOfIndex(index + 1);
            double dyDx = (_data[index + 1] - _data[index]) / (xRight - xLeft);
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
                throw new Exception($"input value b = {b} out of range [{_xmin},{_xmax}]");
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

            int indexOfa = (int)IndexOfX(a);
            int indexOfb = (int)IndexOfX(b);
            if (indexOfa < 0 || indexOfa > _data.Length - 1)
            {
                throw new Exception("left boundary out of range");
            }

            if (indexOfb < 0 || indexOfb > _data.Length - 1)
            {
                throw new Exception("right boundary out of range");
            }

            double result;
            if (indexOfa == indexOfb)
            {
                double xLeft = XOfIndex(indexOfa);
                double xRight = XOfIndex(indexOfa + 1);
                double dyDx = (_data[indexOfa + 1] - _data[indexOfa]) / (xRight - xLeft);
                double ya = _data[indexOfa] + dyDx * (a - xLeft);
                double yb = _data[indexOfa] + dyDx * (b - xLeft);
                result = 0.5 * (ya + yb) * (b - a);
            }
            else
            {
                result = 0.0;
                result += 0.5 * (_data[indexOfa + 1] + Evaluate(a)) * (XOfIndex(indexOfa + 1) - a);
                result += 0.5 * (Evaluate(b) + _data[indexOfb]) * (b - XOfIndex(indexOfb));

                if ((indexOfb - indexOfa) > 1)
                {
                    for (int i = indexOfa + 1; i < indexOfb; i++)
                    {
                        result += 0.5 * (_data[i] + _data[i + 1]) * (XOfIndex(i + 1) - XOfIndex(i));
                    }
                }
            }

            return result * signum;
        }

        /// <summary>
        /// 返回节点索引对应的 x 值。
        /// </summary>
        public double XOfIndex(int index)
        {
            double xs = MathCoreInfo.Pi * index / (_data.Length - 1.0);
            if (_mode > 0)
            {
                xs *= 0.5;
            }

            double fu;
            if (_mode == 2)
            {
                fu = Math.Sin(xs);
            }
            else
            {
                fu = 1.0 - Math.Cos(xs);
                if (_mode == 0)
                {
                    fu *= 0.5;
                }
            }

            return _xmin + (_xmax - _xmin) * fu;
        }

        /// <summary>
        /// 返回给定 x 对应的节点索引。
        /// </summary>
        public double IndexOfX(double x)
        {
            double fu = (x - _xmin) / (_xmax - _xmin);
            fu = Math.Min(1.0, Math.Max(0.0, fu));
            double xs;
            if (_mode == 2)
            {
                xs = Math.Asin(fu);
            }
            else if (_mode == 1)
            {
                xs = Math.Acos(1.0 - fu);
            }
            else
            {
                xs = Math.Acos(1.0 - 2.0 * fu);
            }

            double result = xs / MathCoreInfo.Pi * (_data.Length - 1);
            if (_mode > 0)
            {
                result *= 2.0;
            }

            return result;
        }
    }
}
