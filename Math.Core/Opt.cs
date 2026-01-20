using System;

namespace MathTools.Core
{
    /// <summary>
    /// OPT 优化工具箱基础类（版本号）。
    /// </summary>
    public static class OPT
    {
        /// <summary>
        /// 获取版本号字符串。
        /// </summary>
        public static string GetVersion()
        {
            return "OPT_V03.01.04";
        }
    }

    /// <summary>
    /// Optimization Tool Box (SB1SQPA0) 模块说明（中文翻译）。
    /// </summary>
    /// <remarks>
    /// <para>1. 优化问题分类：</para>
    /// <list type="bullet">
    /// <item>无约束单变量最小化：没有约束，且目标函数为单变量。</item>
    /// <item>线性规划：目标函数与约束为仿射函数。</item>
    /// <item>二次规划：目标函数为二次型，约束为仿射函数。</item>
    /// <item>非线性规划：目标函数与约束均为一般非线性函数。</item>
    /// <item>混合整数非线性规划：部分变量为整数。</item>
    /// </list>
    /// <para>2. 线性规划：</para>
    /// <para>使用修正单纯形法（simplex），并提供 linprog 封装以简化输入；适用于中小规模问题。</para>
    /// <para>3. 二次规划：</para>
    /// <para>使用主动集法 quadprog，并通过 simplex 寻找可行初始点；适用于中小规模问题。</para>
    /// <para>4. 非线性规划：</para>
    /// <para>使用 SQP（Sequential Quadratic Programming），在迭代过程中调用 quadprog，并用 simplex 做可行性初始化；适用于中小规模问题。</para>
    /// <para>5. 混合整数非线性规划：</para>
    /// <para>使用分支定界（BnB）算法，节点处调用 SQP 解决放松问题。请注意该类问题规模扩大后会显著增加时间与内存开销。</para>
    /// <para>6. 测试示例：</para>
    /// <list type="bullet">
    /// <item>单变量最小化：test1Dmin.cpp</item>
    /// <item>线性规划：testLP.cpp</item>
    /// <item>二次规划：testQP.cpp + qpData.txt</item>
    /// <item>非线性规划：testSQP.cpp + exampleNonLinProg.h</item>
    /// <item>混合整数非线性规划：testBnB.cpp</item>
    /// </list>
    /// </remarks>
    public static class OptDocumentation
    {
    }
}
