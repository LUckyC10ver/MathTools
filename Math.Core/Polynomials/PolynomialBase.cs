namespace Math.Core.Polynomials
{
    /// <summary>
    /// 多项式基类。
    /// </summary>
    public abstract class PolynomialBase : DifferentiableFunction
    {
        /// <summary>
        /// 多项式阶数。
        /// </summary>
        public abstract int Degree { get; }

        /// <summary>
        /// 返回导数函数对象。
        /// </summary>
        public override abstract UnaryFunction GetDerivativeFunction();

        /// <summary>
        /// 返回指定点的导数值。
        /// </summary>
        public abstract double GetDerivativeValue(double argument);

        /// <summary>
        /// 计算多项式值。
        /// </summary>
        public override abstract double Invoke(double argument);

        /// <summary>
        /// 设置系数。
        /// </summary>
        public abstract void SetCoefficient(int index, double coeff);
    }
}
