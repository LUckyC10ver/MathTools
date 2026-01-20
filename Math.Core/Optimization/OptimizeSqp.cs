using System;
using System.IO;

namespace MathTools.Core
{
    public class SqpOptions
    {
        public int Display { get; set; }
        public double TolArg { get; set; }
        public double TolObj { get; set; }
        public double TolCon { get; set; }
        public int LineSearch { get; set; }
        public bool CheckGrad { get; set; }
        public int MaxFunEvals { get; set; }
        public double DiffMinChange { get; set; }
        public double DiffMaxChange { get; set; }
        public double MaxDouble { get; set; }
        public double Eps { get; set; }
        public TextWriter SqpLog { get; set; }
        public TextWriter QpLog { get; set; }
        public int MaxNumItersQp { get; set; }

        public SqpOptions(
            int display = 0,
            TextWriter sqpLog = null,
            double tolArg = 1e-2,
            double tolObj = 1e-4,
            double tolCon = 1e-6,
            double eps = 1e-11,
            bool checkGrad = false,
            int maxFunEvals = 0,
            double diffMinChange = 1e-9,
            double diffMaxChange = 1e-6,
            double maxDouble = 1.7e38,
            TextWriter qpLog = null,
            int lineSearch = 0,
            int maxItersQp = 1000)
        {
            if (lineSearch < 0 || lineSearch > 5)
            {
                throw new Exception($"wrong line search option: {lineSearch}; must be >= 0 and <= 5.");
            }

            Display = display;
            TolArg = tolArg;
            TolObj = tolObj;
            TolCon = tolCon;
            LineSearch = lineSearch;
            CheckGrad = checkGrad;
            MaxFunEvals = maxFunEvals;
            DiffMinChange = diffMinChange;
            DiffMaxChange = diffMaxChange;
            MaxDouble = maxDouble;
            Eps = eps;
            SqpLog = sqpLog;
            QpLog = qpLog;
            MaxNumItersQp = maxItersQp;
        }
    }

    public class SqpInfo
    {
        public double ObjValue { get; set; }
        public int SqpCount { get; set; }
        public int FunCount { get; set; }
        public int GradCount { get; set; }
        public double StepLength { get; set; }
        public string HowQp { get; set; }
        public string How { get; set; }

        public SqpInfo()
        {
            ObjValue = 0.0;
            SqpCount = 0;
            FunCount = 0;
            GradCount = 0;
            StepLength = 1.0;
            HowQp = string.Empty;
            How = string.Empty;
        }
    }

    public partial class Functions
    {
        public delegate void SqpObjective(out double f, out double[] grad, double[] x, bool calcGrad);
        public delegate void SqpConstraints(out int nequ, out double[] g, out double[][] grad, double[] x, bool calcGrad);

        public static void runSqp(
            out double[] xout,
            out SqpInfo info,
            out double[] lambda,
            out int[] actSet,
            SqpObjective obj,
            double[] lb,
            double[] ub,
            double[][] Aequ,
            double[] bequ,
            double[][] Aieq,
            double[] bieq,
            double[] xOrig,
            SqpConstraints con,
            SqpOptions opt = null,
            double[][] Hess = null)
        {
            if (obj == null)
            {
                throw new Exception("runSqp: objective is required");
            }

            if (con == null)
            {
                throw new Exception("runSqp: constraints are required");
            }

            if (xOrig == null || xOrig.Length == 0)
            {
                throw new Exception("runSqp: starting point is required");
            }

            opt ??= new SqpOptions();
            info = new SqpInfo();
            xout = (double[])xOrig.Clone();
            lambda = Array.Empty<double>();
            actSet = Array.Empty<int>();

            int nvars = xout.Length;
            double[] lower = lb ?? Array.Empty<double>();
            double[] upper = ub ?? Array.Empty<double>();

            if (lower.Length > 0 && lower.Length != nvars)
            {
                throw new Exception("runSqp: lb size mismatch");
            }

            if (upper.Length > 0 && upper.Length != nvars)
            {
                throw new Exception("runSqp: ub size mismatch");
            }

            if ((Aequ == null || Aequ.Length == 0) && bequ != null && bequ.Length > 0)
            {
                throw new Exception("runSqp: Aequ and bequ are incompatible");
            }

            if ((Aieq == null || Aieq.Length == 0) && bieq != null && bieq.Length > 0)
            {
                throw new Exception("runSqp: Aieq and bieq are incompatible");
            }

            if (Aequ == null)
            {
                Aequ = Array.Empty<double[]>();
            }

            if (Aieq == null)
            {
                Aieq = Array.Empty<double[]>();
            }

            bequ ??= Array.Empty<double>();
            bieq ??= Array.Empty<double>();

            if (Aequ.Length > 0 && Aequ[0].Length != nvars)
            {
                throw new Exception("runSqp: Aequ column size mismatch");
            }

            if (Aieq.Length > 0 && Aieq[0].Length != nvars)
            {
                throw new Exception("runSqp: Aieq column size mismatch");
            }

            int lenvlb = lower.Length;
            int lenvub = upper.Length;
            for (int i = 0; i < Math.Min(lenvlb, lenvub); i++)
            {
                if (lower[i] > upper[i])
                {
                    throw new Exception("Bounds Infeasible");
                }
            }

            for (int i = 0; i < lenvlb; i++)
            {
                if (xout[i] < lower[i])
                {
                    xout[i] = lower[i] + opt.Eps;
                }
            }

            for (int i = 0; i < lenvub; i++)
            {
                if (xout[i] > upper[i])
                {
                    xout[i] = upper[i];
                }
            }

            if (Hess == null || Hess.Length != nvars || Hess[0].Length != nvars)
            {
                Hess = CreateMatrix(nvars, nvars, 0.0);
                unit_matrix(ref Hess);
            }

            if (opt.Display > 0 && opt.SqpLog != null)
            {
                opt.SqpLog.WriteLine("version number of runSqp: 1.0");
            }

            obj(out double f, out double[] gf, xout, true);
            bool gftype = gf != null && gf.Length > 0;

            con(out int nequ, out double[] g, out double[][] AN, xout, true);
            int ncstr = g?.Length ?? 0;
            bool ggtype = AN != null && AN.Length > 0;

            if (g == null)
            {
                g = Array.Empty<double>();
            }

            if (AN == null)
            {
                AN = Array.Empty<double[]>();
            }

            if (ncstr == 0)
            {
                nequ = 0;
                AN = CreateMatrix(0, nvars, 0.0);
            }

            if (opt.MaxFunEvals == 0)
            {
                int baseCount = (int)Math.Log(nvars + Math.Exp(1.0)) * 300;
                opt.MaxFunEvals = gftype && ggtype ? baseCount : baseCount * nvars;
            }

            info.FunCount = 1;
            info.GradCount = 1;
            info.StepLength = 1.0;

            double bestf = opt.MaxDouble;
            double[] bestx = (double[])xout.Clone();

            double[] oldx = (double[])xout.Clone();
            double[] oldg = (double[])g.Clone();
            double[] oldgf = gftype ? (double[])gf.Clone() : new double[nvars];
            double[][] oldAN = AN.Length == 0 ? CreateMatrix(0, nvars, 0.0) : CloneMatrix(AN);

            double[] lambdaNonlinear = new double[ncstr];
            double[] oldLambdaNonlinear = new double[ncstr];
            double[] lambdaLinear = new double[bequ.Length + bieq.Length];
            double[] oldLambdaLinear = new double[lambdaLinear.Length];

            if (opt.Display > 0 && opt.SqpLog != null)
            {
                opt.SqpLog.WriteLine("\n iter   f-COUNT   FUNCTION       MAX{g}         STEP  Procedures");
            }

            int status = 0;
            while (status != 1 && info.SqpCount < opt.MaxFunEvals)
            {
                info.SqpCount++;

                if (!gftype || (!ggtype && ncstr > 0) || opt.CheckGrad)
                {
                    double origf = f;
                    double[] origg = (double[])g.Clone();
                    int diffN = nvars;
                    double[] gfFD = new double[nvars];
                    double[][] ggFD = CreateMatrix(ncstr, nvars, 0.0);

                    for (int gcnt = 0; gcnt < diffN; gcnt++)
                    {
                        double dx = Math.Abs(xout[gcnt]) * 1e-7 + 1e-7;
                        dx = Math.Max(dx, opt.DiffMinChange);
                        dx = Math.Min(dx, opt.DiffMaxChange);

                        if (lenvlb > 0 && lenvub > 0)
                        {
                            if (xout[gcnt] + dx > upper[gcnt] && xout[gcnt] - dx >= lower[gcnt])
                            {
                                dx = -dx;
                            }
                            else if (xout[gcnt] + dx > upper[gcnt] && xout[gcnt] - dx < lower[gcnt])
                            {
                                double dxMaxO = upper[gcnt] - xout[gcnt];
                                double dxMaxU = xout[gcnt] - lower[gcnt];
                                if (dxMaxO > dxMaxU)
                                {
                                    dx = dxMaxO;
                                }
                                else if (dxMaxU > opt.Eps)
                                {
                                    dx = -dxMaxU;
                                }
                                else
                                {
                                    continue;
                                }
                            }
                        }
                        else if (lenvlb > 0 && xout[gcnt] + dx < lower[gcnt])
                        {
                            if (xout[gcnt] - dx >= lower[gcnt])
                            {
                                dx = -dx;
                            }
                            else
                            {
                                throw new Exception($"{gcnt + 1}-th variable: lower bound > initial value");
                            }
                        }
                        else if (lenvub > 0 && xout[gcnt] + dx > upper[gcnt])
                        {
                            if (xout[gcnt] - dx <= upper[gcnt])
                            {
                                dx = -dx;
                            }
                            else
                            {
                                throw new Exception($"{gcnt + 1}-th variable: upper bound < initial value");
                            }
                        }

                        double xo = xout[gcnt];
                        xout[gcnt] += dx;

                        if (!gftype || opt.CheckGrad)
                        {
                            obj(out f, out _, xout, false);
                        }

                        if (ncstr > 0 && (!ggtype || opt.CheckGrad))
                        {
                            con(out nequ, out g, out _, xout, false);
                        }

                        gfFD[gcnt] = (f - origf) / dx;
                        if (ncstr > 0 && (!ggtype || opt.CheckGrad))
                        {
                            int diff = g.Length;
                            for (int i = 0; i < diff; i++)
                            {
                                ggFD[i][gcnt] = (g[i] - origg[i]) / dx;
                            }
                        }

                        xout[gcnt] = xo;
                    }

                    if (!gftype || opt.CheckGrad)
                    {
                        gf = gfFD;
                        gftype = true;
                    }

                    if (ncstr > 0 && (!ggtype || opt.CheckGrad))
                    {
                        AN = ggFD;
                        ggtype = true;
                    }

                    info.FunCount += nvars;
                    f = origf;
                    g = origg;
                }

                if (info.FunCount > 1)
                {
                    if (gftype)
                    {
                        obj(out f, out gf, xout, true);
                    }

                    if (ggtype && ncstr > 0)
                    {
                        con(out nequ, out g, out AN, xout, true);
                    }
                }

                info.How = string.Empty;

                for (int i = 0; i < nequ; i++)
                {
                    if (dot(AN[i], gf) > 0.0)
                    {
                        for (int j = 0; j < nvars; j++)
                        {
                            AN[i][j] *= -1.0;
                        }

                        g[i] *= -1.0;
                    }
                }

                if (info.GradCount > 1)
                {
                    double[] gnew = new double[nvars];
                    if (lambdaNonlinear.Length == 0)
                    {
                        Array.Copy(gf, gnew, nvars);
                    }
                    else
                    {
                        double[] tmp = vecdotmat(AN, lambdaNonlinear);
                        for (int i = 0; i < nvars; i++)
                        {
                            gnew[i] = gf[i] + tmp[i];
                        }
                    }

                    double[] gold = new double[nvars];
                    if (oldAN.Length > 0)
                    {
                        double[] tmp = vecdotmat(oldAN, oldLambdaNonlinear);
                        for (int i = 0; i < nvars; i++)
                        {
                            gold[i] = oldgf[i] + tmp[i];
                        }
                    }
                    else
                    {
                        Array.Copy(oldgf, gold, nvars);
                    }

                    double[] q = new double[nvars];
                    double[] s = new double[nvars];
                    for (int i = 0; i < nvars; i++)
                    {
                        q[i] = gnew[i] - gold[i];
                        s[i] = xout[i] - oldx[i];
                    }

                    if (dot(q, s) < info.StepLength * info.StepLength * 1e-3)
                    {
                        int count = 0;
                        while (dot(q, s) < -1e-5 && count++ < 1000)
                        {
                            int yind = nvars - 1;
                            double pMin = q[yind] * s[yind];
                            for (int i = nvars - 2; i >= 0; i--)
                            {
                                double p = q[i] * s[i];
                                if (p < pMin)
                                {
                                    yind = i;
                                    pMin = p;
                                }
                            }

                            q[yind] /= 2.0;
                        }

                        if (dot(q, s) < opt.Eps * normFrobenius(Hess))
                        {
                            double[] gAN = new double[nvars];
                            if (g.Length == 1)
                            {
                                for (int i = 0; i < nvars; i++)
                                {
                                    gAN[i] = AN[0][i] * g[0];
                                }
                            }
                            else if (g.Length > 1)
                            {
                                gAN = vecdotmat(AN, g);
                            }

                            double[] gANold = new double[nvars];
                            if (oldg.Length == 1)
                            {
                                for (int i = 0; i < nvars; i++)
                                {
                                    gANold[i] = oldAN[0][i] * oldg[0];
                                }
                            }
                            else if (oldg.Length > 1)
                            {
                                gANold = vecdotmat(oldAN, oldg);
                            }

                            double[] factor = new double[nvars];
                            double maxv = 0.0;
                            for (int i = 0; i < nvars; i++)
                            {
                                double dgAN = gAN[i] - gANold[i];
                                if (s[i] * dgAN > 0.0 && s[i] * q[i] <= opt.Eps)
                                {
                                    factor[i] = dgAN;
                                    double absVal = Math.Abs(dgAN);
                                    if (absVal > maxv)
                                    {
                                        maxv = absVal;
                                    }
                                }
                            }

                            if (maxv < opt.Eps)
                            {
                                for (int i = 0; i < nvars; i++)
                                {
                                    factor[i] = 1e-5 * Sign(s[i]);
                                }
                            }

                            double wt = 1e-2;
                            while (dot(q, s) < opt.Eps * normFrobenius(Hess) && wt < 1.0 / opt.Eps)
                            {
                                for (int i = 0; i < nvars; i++)
                                {
                                    q[i] += factor[i] * wt;
                                }

                                wt *= 2.0;
                            }
                        }
                    }

                    double qs = dot(q, s);
                    if (qs > opt.Eps)
                    {
                        double[] Hs = matdotvec(Hess, s);
                        double sHs = dot(s, Hs);
                        if (sHs <= 1.0 / opt.MaxDouble)
                        {
                            throw new Exception($"error in BFGS update: last hessian H is not positive definite, s'Hs= {sHs}");
                        }

                        double alpha = 1.0 / qs;
                        double beta = -1.0 / sHs;
                        rank1Update(ref Hess, q, alpha);
                        rank1Update(ref Hess, Hs, beta);
                    }
                }
                else
                {
                    double ngf = opt.Eps + norm(gf);
                    for (int i = 0; i < ncstr; i++)
                    {
                        oldLambdaNonlinear[i] = ngf / (norm(AN[i]) + opt.Eps);
                    }
                    for (int i = 0; i < bequ.Length; i++)
                    {
                        oldLambdaLinear[i] = ngf / (norm(Aequ[i]) + opt.Eps);
                    }
                    for (int i = 0; i < bieq.Length; i++)
                    {
                        oldLambdaLinear[i + bequ.Length] = ngf / (norm(Aieq[i]) + opt.Eps);
                    }
                }

                info.GradCount++;
                oldLambdaNonlinear = (double[])lambdaNonlinear.Clone();
                oldAN = CloneMatrix(AN);
                oldgf = (double[])gf.Clone();
                oldg = (double[])g.Clone();
                double oldf = f;
                oldx = (double[])xout.Clone();

                double[] sd = new double[nvars];
                int nLinEqu = bequ.Length;
                int nLinIeq = bieq.Length;
                int qpNequ = nequ + nLinEqu;
                int nieq = ncstr - nequ;
                int qpNieq = nieq + nLinIeq;

                double[][] qpAequ = CreateMatrix(qpNequ, nvars, 0.0);
                double[] qpBequ = new double[qpNequ];
                double[][] qpAieq = CreateMatrix(qpNieq, nvars, 0.0);
                double[] qpBieq = new double[qpNieq];
                double[] qpLb = new double[lenvlb];
                double[] qpUb = new double[lenvub];

                for (int i = 0; i < lenvlb; i++)
                {
                    qpLb[i] = lower[i] - xout[i];
                }
                for (int i = 0; i < lenvub; i++)
                {
                    qpUb[i] = upper[i] - xout[i];
                }

                if (nLinEqu > 0)
                {
                    double[] tmpEq = matdotvec(Aequ, xout);
                    for (int i = 0; i < nLinEqu; i++)
                    {
                        qpAequ[i] = (double[])Aequ[i].Clone();
                        qpBequ[i] = bequ[i] - tmpEq[i];
                    }
                }

                for (int i = 0; i < nequ; i++)
                {
                    qpAequ[i + nLinEqu] = (double[])AN[i].Clone();
                    qpBequ[i + nLinEqu] = -g[i];
                }

                for (int i = 0; i < nieq; i++)
                {
                    qpAieq[i] = (double[])AN[i + nequ].Clone();
                    qpBieq[i] = -g[i + nequ];
                }

                if (nLinIeq > 0)
                {
                    double[] tmpIeq = matdotvec(Aieq, xout);
                    for (int i = 0; i < nLinIeq; i++)
                    {
                        qpAieq[i + nieq] = (double[])Aieq[i].Clone();
                        qpBieq[i + nieq] = bieq[i] - tmpIeq[i];
                    }
                }

                int qpErr = 0;
                info.HowQp = string.Empty;
                try
                {
                    int verb = opt.QpLog != null && opt.Display > 3 ? opt.Display - 3 : 0;
                    quadprog(out sd, out double[] lambdaQp, out int[] actQp, out _, Hess, gf, qpAieq, qpBieq, qpAequ, qpBequ, qpLb, qpUb, new double[nvars], verb, opt.Eps, opt.MaxDouble, opt.MaxNumItersQp);

                    lambda = MapQpLambda(lambdaQp, nequ, nLinEqu, nLinIeq, nieq, lenvlb, lenvub);
                    actSet = actQp;
                    lambdaNonlinear = new double[ncstr];
                    for (int i = 0; i < ncstr; i++)
                    {
                        int offset = nLinEqu + i;
                        if (offset < lambdaQp.Length)
                        {
                            lambdaNonlinear[i] = lambdaQp[offset];
                        }
                    }

                    lambdaLinear = new double[nLinEqu + nLinIeq];
                    for (int i = 0; i < nLinEqu && i < lambdaQp.Length; i++)
                    {
                        lambdaLinear[i] = lambdaQp[i];
                    }
                    for (int i = 0; i < nLinIeq; i++)
                    {
                        int offset = nLinEqu + nequ + nieq + i;
                        if (offset < lambdaQp.Length)
                        {
                            lambdaLinear[nLinEqu + i] = lambdaQp[offset];
                        }
                    }
                }
                catch (Exception ex)
                {
                    qpErr = 1;
                    info.HowQp = ex.Message;
                    lambda = new double[nequ + nLinEqu + nLinIeq + nieq + lenvlb + lenvub];
                    actSet = Array.Empty<int>();
                }

                double[] ga = new double[ncstr];
                double mg = -opt.MaxDouble;
                for (int i = 0; i < nequ; i++)
                {
                    double d = Math.Abs(g[i]);
                    ga[i] = d;
                    if (d > mg)
                    {
                        mg = d;
                    }
                }
                for (int i = 0; i < nieq; i++)
                {
                    double d = g[i + nequ];
                    ga[i + nequ] = d;
                    if (d > mg)
                    {
                        mg = d;
                    }
                }

                double mgLin = -opt.MaxDouble;
                double[] gaLin = new double[nLinEqu + nLinIeq];
                for (int i = 0; i < nLinEqu; i++)
                {
                    double d = Math.Abs(qpBequ[i]);
                    gaLin[i] = d;
                    if (d > mgLin)
                    {
                        mgLin = d;
                    }
                }
                for (int i = 0; i < nLinIeq; i++)
                {
                    double d = -qpBieq[i + nieq];
                    gaLin[i + nLinEqu] = d;
                    if (d > mgLin)
                    {
                        mgLin = d;
                    }
                }

                mg = Math.Max(mg, mgLin);

                for (int i = 0; i < ncstr; i++)
                {
                    double mid = 0.5 * (lambdaNonlinear[i] + oldLambdaNonlinear[i]);
                    oldLambdaNonlinear[i] = Math.Max(mid, lambdaNonlinear[i]);
                }
                for (int i = 0; i < nLinEqu + nLinIeq; i++)
                {
                    double mid = 0.5 * (lambdaLinear[i] + oldLambdaLinear[i]);
                    oldLambdaLinear[i] = Math.Max(mid, lambdaLinear[i]);
                }

                double[] matx = (double[])xout.Clone();
                double MATL = Merit(f, ga, oldLambdaNonlinear, gaLin, oldLambdaLinear);
                double MERIT = opt.MaxDouble;
                info.StepLength = 2.0;

                while (MERIT > MATL && info.FunCount < opt.MaxFunEvals)
                {
                    info.StepLength /= 2.0;
                    if (Math.Abs(info.StepLength) < 1e-4)
                    {
                        info.StepLength *= -1.0;
                    }

                    for (int i = 0; i < nvars; i++)
                    {
                        xout[i] = matx[i] + sd[i] * info.StepLength;
                    }

                    obj(out f, out _, xout, false);
                    if (ncstr > 0)
                    {
                        con(out nequ, out g, out _, xout, false);
                    }
                    info.FunCount++;

                    mg = -opt.MaxDouble;
                    for (int i = 0; i < nequ; i++)
                    {
                        double d = Math.Abs(g[i]);
                        ga[i] = d;
                        if (d > mg)
                        {
                            mg = d;
                        }
                    }
                    for (int i = 0; i < nieq; i++)
                    {
                        double d = g[i + nequ];
                        ga[i + nequ] = d;
                        if (d > mg)
                        {
                            mg = d;
                        }
                    }

                    mgLin = -opt.MaxDouble;
                    for (int i = 0; i < nLinEqu; i++)
                    {
                        double d = Math.Abs(dot(Aequ[i], xout) - bequ[i]);
                        gaLin[i] = d;
                        if (d > mgLin)
                        {
                            mgLin = d;
                        }
                    }
                    for (int i = 0; i < nLinIeq; i++)
                    {
                        double d = dot(Aieq[i], xout) - bieq[i];
                        gaLin[i + nLinEqu] = d;
                        if (d > mgLin)
                        {
                            mgLin = d;
                        }
                    }

                    mg = Math.Max(mg, mgLin);
                    MERIT = Merit(f, ga, oldLambdaNonlinear, gaLin, oldLambdaLinear);
                }

                if (mg < opt.TolCon && f < bestf)
                {
                    bestf = f;
                    bestx = (double[])xout.Clone();
                }

                double help1 = 0.0;
                for (int i = 0; i < sd.Length; i++)
                {
                    help1 = Math.Max(Math.Abs(sd[i]), help1);
                }

                double help2 = Math.Abs(dot(gf, sd));

                if (opt.Display > 0 && opt.SqpLog != null)
                {
                    string qpmsg = qpErr == 0 ? " " : info.HowQp;
                    opt.SqpLog.WriteLine($"\n {info.SqpCount,5} {info.FunCount,5} {f,12} {mg,12} {info.StepLength,12} {info.How}    {qpmsg}");
                }

                if (help1 < 2.0 * opt.TolArg && help2 < 2.0 * opt.TolObj && mg < opt.TolCon)
                {
                    if (bestf < f)
                    {
                        xout = (double[])bestx.Clone();
                        info.ObjValue = bestf;
                    }
                    else
                    {
                        info.ObjValue = f;
                    }

                    status = 1;
                }
                else if (info.FunCount >= opt.MaxFunEvals)
                {
                    if (oldf > bestf)
                    {
                        xout = (double[])bestx.Clone();
                        info.ObjValue = bestf;
                    }
                    else
                    {
                        xout = (double[])matx.Clone();
                        info.ObjValue = oldf;
                    }

                    throw new Exception("Maximum number of iterations exceeded");
                }
            }

            info.ObjValue = f;
            if (bestf < f)
            {
                xout = (double[])bestx.Clone();
                info.ObjValue = bestf;
            }

            if (lambda == null)
            {
                lambda = Array.Empty<double>();
            }

            if (actSet == null)
            {
                actSet = Array.Empty<int>();
            }
        }

        private static double[][] CreateMatrix(int rows, int cols, double value)
        {
            return new double[0][].Resize(rows, cols, value);
        }

        private static double[][] CloneMatrix(double[][] src)
        {
            var result = CreateMatrix(src.Length, src.Length == 0 ? 0 : src[0].Length, 0.0);
            for (int i = 0; i < src.Length; i++)
            {
                Array.Copy(src[i], result[i], src[i].Length);
            }

            return result;
        }

        private static double Sign(double value)
        {
            return value >= 0 ? 1.0 : -1.0;
        }

        private static double Merit(
            double f,
            double[] ga,
            double[] lambda,
            double[] gaLin,
            double[] lambdaLin)
        {
            double merit = f + 1e-30;
            for (int i = 0; i < ga.Length; i++)
            {
                merit += Math.Max(ga[i], 0.0) * lambda[i];
            }
            for (int i = 0; i < gaLin.Length; i++)
            {
                merit += Math.Max(gaLin[i], 0.0) * lambdaLin[i];
            }

            return merit;
        }

        private static double[] MapQpLambda(
            double[] lambdaQp,
            int nNonlinEqu,
            int nLinEqu,
            int nLinIeq,
            int nNonlinIeq,
            int nLb,
            int nUb)
        {
            int total = nNonlinEqu + nLinEqu + nLinIeq + nNonlinIeq + nLb + nUb;
            double[] lambda = new double[total];

            int qpEq = nLinEqu + nNonlinEqu;
            int qpIeq = nNonlinIeq + nLinIeq;
            int qpOffsetIeq = qpEq;

            for (int i = 0; i < nNonlinEqu; i++)
            {
                int src = nLinEqu + i;
                if (src < lambdaQp.Length)
                {
                    lambda[i] = Math.Abs(lambdaQp[src]);
                }
            }

            for (int i = 0; i < nLinEqu; i++)
            {
                if (i < lambdaQp.Length)
                {
                    lambda[nNonlinEqu + i] = Math.Abs(lambdaQp[i]);
                }
            }

            for (int i = 0; i < nLinIeq; i++)
            {
                int src = qpOffsetIeq + nNonlinIeq + i;
                int dst = nNonlinEqu + nLinEqu + i;
                if (src < lambdaQp.Length)
                {
                    lambda[dst] = lambdaQp[src];
                }
            }

            for (int i = 0; i < nNonlinIeq; i++)
            {
                int src = qpOffsetIeq + i;
                int dst = nNonlinEqu + nLinEqu + nLinIeq + i;
                if (src < lambdaQp.Length)
                {
                    lambda[dst] = lambdaQp[src];
                }
            }

            int baseIndex = nNonlinEqu + nLinEqu + nLinIeq + nNonlinIeq;
            int qpBoundStart = qpEq + qpIeq;
            for (int i = 0; i < nLb; i++)
            {
                int src = qpBoundStart + i;
                if (src < lambdaQp.Length)
                {
                    lambda[baseIndex + i] = lambdaQp[src];
                }
            }

            for (int i = 0; i < nUb; i++)
            {
                int src = qpBoundStart + nLb + i;
                if (src < lambdaQp.Length)
                {
                    lambda[baseIndex + nLb + i] = lambdaQp[src];
                }
            }

            return lambda;
        }
    }
}
