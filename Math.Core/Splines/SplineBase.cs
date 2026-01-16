using System;

namespace MathTools.Core.Splines
{
    /// <summary>
    /// 样条基类，负责控制点管理与排序。
    /// </summary>
    public abstract class SplineBase
    {
        protected double[] _xValues = Array.Empty<double>();
        protected double[] _yValues = Array.Empty<double>();
        protected bool _changed = true;
        protected double _epsilon = MathCoreInfo.BcEpsEquidistantLinearSpline;

        /// <summary>
        /// 控制点数量。
        /// </summary>
        public int Size => _xValues.Length;

        /// <summary>
        /// 控制点 X 值。
        /// </summary>
        public double[] XValues => _xValues;

        /// <summary>
        /// 控制点 Y 值。
        /// </summary>
        public double[] YValues => _yValues;

        /// <summary>
        /// 初始化控制点，值填充为固定值。
        /// </summary>
        public void Init(double[] points, double initValue = 0.0, bool sorted = false)
        {
            _xValues = points ?? throw new Exception("points is null");
            _yValues = new double[_xValues.Length];
            for (int i = 0; i < _yValues.Length; i++)
            {
                _yValues[i] = initValue;
            }

            if (!sorted)
            {
                Functions.sortPoints(ref _xValues, ref _yValues);
            }

            _changed = true;
            UpdateEpsilon();
        }

        /// <summary>
        /// 初始化控制点与取值。
        /// </summary>
        public void Init(double[] points, double[] values, bool sorted = false)
        {
            if (points == null || values == null || points.Length != values.Length)
            {
                throw new Exception("init: vector sizes differ");
            }

            _xValues = (double[])points.Clone();
            _yValues = (double[])values.Clone();

            if (!sorted)
            {
                Functions.sortPoints(ref _xValues, ref _yValues);
            }

            _changed = true;
            UpdateEpsilon();
        }

        /// <summary>
        /// 设置值数组。
        /// </summary>
        public void SetValues(double[] values)
        {
            if (values == null || values.Length != _yValues.Length)
            {
                throw new Exception("setValues: vector sizes differ");
            }

            _yValues = (double[])values.Clone();
            _changed = true;
        }

        /// <summary>
        /// X 最小值。
        /// </summary>
        public double XMin()
        {
            if (_xValues.Length == 0)
            {
                throw new Exception("Range not set");
            }

            return _xValues[0];
        }

        /// <summary>
        /// X 最大值。
        /// </summary>
        public double XMax()
        {
            if (_xValues.Length == 0)
            {
                throw new Exception("Range not set");
            }

            return _xValues[_xValues.Length - 1];
        }

        /// <summary>
        /// 是否有修改。
        /// </summary>
        public bool IsChanged() => _changed;

        /// <summary>
        /// 重置修改标记。
        /// </summary>
        public void ResetChanged()
        {
            _changed = false;
        }

        /// <summary>
        /// 更新内部数据。
        /// </summary>
        protected virtual void Update()
        {
            if (_changed)
            {
                _changed = false;
            }
        }

        /// <summary>
        /// 查找小于等于 x 的最大索引。
        /// </summary>
        protected int IndexFromValue(double x)
        {
            int n = _xValues.Length;
            int klo = 0;
            int khi = n;
            while (khi - klo > 1)
            {
                int k = (khi + klo) / 2;
                if (_xValues[k] > x)
                {
                    khi = k;
                }
                else
                {
                    klo = k;
                }
            }

            return klo;
        }

        private void UpdateEpsilon()
        {
            if (_xValues.Length == 0)
            {
                _epsilon = MathCoreInfo.BcEpsEquidistantLinearSpline;
                return;
            }

            double maxAbs = Math.Max(Math.Abs(_xValues[0]), Math.Abs(_xValues[^1]));
            _epsilon = MathCoreInfo.BcEpsEquidistantLinearSpline * maxAbs;
        }
    }
}
