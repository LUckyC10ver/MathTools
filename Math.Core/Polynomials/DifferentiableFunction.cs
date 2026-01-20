namespace MathTools.Core
{
    /// <summary>
    /// 可微一元函数基类。
    /// </summary>
    public abstract class DifferentiableFunction : UnaryFunction
    {
        /// <summary>
        /// 返回导数函数对象（默认无）。
        /// </summary>
        public virtual UnaryFunction GetDerivativeFunction()
        {
            return null;
        }
    }
}
