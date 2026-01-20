using System;
using System.Collections.Generic;

namespace MathTools.Core
{
    public partial class Functions
    {
        public static int quadprog(
            out double[] x,
            out double[] lambda,
            out int[] actSet,
            out double[][] tableau,
            double[][] H,
            double[] f,
            double[][] Aieq,
            double[] bieq,
            double[][] Aeq,
            double[] beq,
            double[] lbd,
            double[] ubd,
            double[] x0,
            int verbosity = 0,
            double EPS = 1e-16,
            double BIG = 1e99,
            int MAXITER = 1000)
        {
            if (H == null || f == null)
            {
                throw new Exception("quadprog: H and f are required");
            }

            int nvar = f.Length;
            if (nvar == 0)
            {
                throw new Exception("quadprog: objective has no variables");
            }

            if (H.Length != nvar || H[0].Length != nvar)
            {
                throw new Exception("quadprog: H must be a square matrix of size nvar");
            }

            int nieq = bieq?.Length ?? 0;
            int nequ = beq?.Length ?? 0;

            if (nieq > 0 && Aieq == null)
            {
                throw new Exception("quadprog: Aieq is required for bieq");
            }

            if (nequ > 0 && Aeq == null)
            {
                throw new Exception("quadprog: Aeq is required for beq");
            }

            if (Aieq != null && Aieq.Length > 0 && Aieq[0].Length != nvar)
            {
                throw new Exception("quadprog: Aieq column size must match nvar");
            }

            if (Aeq != null && Aeq.Length > 0 && Aeq[0].Length != nvar)
            {
                throw new Exception("quadprog: Aeq column size must match nvar");
            }

            double[] lower = lbd ?? Array.Empty<double>();
            double[] upper = ubd ?? Array.Empty<double>();
            int lbd_cnt = lower.Length;
            int ubd_cnt = upper.Length;
            if (lbd_cnt > 0 && lbd_cnt != nvar)
            {
                throw new Exception("quadprog: lbd size must match nvar");
            }

            if (ubd_cnt > 0 && ubd_cnt != nvar)
            {
                throw new Exception("quadprog: ubd size must match nvar");
            }

            if (x0 == null || x0.Length != nvar)
            {
                throw new Exception("quadprog: size(f) != size(x0)");
            }

            x = (double[])x0.Clone();

            const double MINDOUBLE = 1e-160;
            double tightToleranceRel = EPS / 100.0;
            double tightToleranceAbs = EPS / 100.0;
            double TOL_REL = EPS;
            double TOL_ABS = EPS;

            var fixInd = new List<int>();
            var lbdInd = new List<int>();
            var ubdInd = new List<int>();
            var infeasConstraint = new List<int>();

            if (lbd_cnt > 0 && ubd_cnt > 0)
            {
                int mn = Math.Min(lbd_cnt, ubd_cnt);
                for (int i = 0; i < mn; i++)
                {
                    if (lower[i] > upper[i])
                    {
                        int offset = nequ + nieq;
                        infeasConstraint.Add(i + offset);
                        infeasConstraint.Add(i + offset + lbd_cnt);
                    }
                    else if (upper[i] - lower[i] <= TOL_REL * Math.Abs(upper[i]))
                    {
                        if (verbosity > 0)
                        {
                            Console.WriteLine($"WARNING: lbd[{i}] == ubd[{i}] == {upper[i]} (converted to equality)");
                        }

                        fixInd.Add(i);
                    }
                }

                if (infeasConstraint.Count > 0)
                {
                    throw new Exception($"quadprog: excluding bounds; stopping the optimization. Indices: {string.Join(",", infeasConstraint)}");
                }
            }

            if (lbd_cnt > 0 || ubd_cnt > 0)
            {
                if (lbd_cnt > 0)
                {
                    lbdInd.Clear();
                    for (int i = 0; i < nvar; i++)
                    {
                        if (lower[i] > (-BIG) && upper.Length > i && upper[i] > lower[i])
                        {
                            lbdInd.Add(i);
                        }
                    }
                }

                if (ubd_cnt > 0)
                {
                    ubdInd.Clear();
                    for (int i = 0; i < nvar; i++)
                    {
                        if (upper[i] < BIG && upper[i] > lower[i])
                        {
                            ubdInd.Add(i);
                        }
                    }
                }
            }

            int nold = nequ + nieq;
            int neqold = nequ;
            int noFixed = fixInd.Count;
            nequ += noFixed;
            int ncst = nequ + nieq + lbdInd.Count + ubdInd.Count;

            double[][] A = CreateMatrix(ncst, nvar, 0.0);
            double[] c = new double[ncst];

            for (int i = 0; i < neqold; i++)
            {
                c[i] = beq[i];
                for (int j = 0; j < nvar; j++)
                {
                    A[i][j] = Aeq[i][j];
                }
            }

            for (int i = 0; i < noFixed; i++)
            {
                int ii = fixInd[i];
                c[neqold + i] = lower[ii];
                for (int j = 0; j < nvar; j++)
                {
                    A[neqold + i][j] = 0.0;
                }
                A[neqold + i][ii] = 1.0;
            }

            for (int i = 0; i < nieq; i++)
            {
                c[i + nequ] = bieq[i];
                for (int j = 0; j < nvar; j++)
                {
                    A[i + nequ][j] = Aieq[i][j];
                }
            }

            int insertIdx = nold + noFixed;
            if (lbdInd.Count > 0)
            {
                for (int j = 0; j < lbdInd.Count; j++)
                {
                    int idx = lbdInd[j];
                    for (int l = 0; l < nvar; l++)
                    {
                        A[insertIdx][l] = 0.0;
                    }
                    A[insertIdx][idx] = -1.0;
                    c[insertIdx] = -lower[idx];
                    insertIdx++;
                }
            }

            if (ubdInd.Count > 0)
            {
                for (int j = 0; j < ubdInd.Count; j++)
                {
                    int idx = ubdInd[j];
                    for (int l = 0; l < nvar; l++)
                    {
                        A[insertIdx][l] = 0.0;
                    }
                    A[insertIdx][idx] = 1.0;
                    c[insertIdx] = upper[idx];
                    insertIdx++;
                }
            }

            double[] normA = new double[ncst];
            for (int i = 0; i < ncst; i++)
            {
                double normVal = 0.0;
                for (int j = 0; j < nvar; j++)
                {
                    double d = A[i][j];
                    normVal += d * d;
                }
                normVal = Math.Sqrt(normVal);
                if (normVal <= TOL_ABS && Math.Abs(c[i]) > TOL_ABS)
                {
                    if (i < nequ)
                    {
                        infeasConstraint.Add(i);
                    }
                    else if (c[i] < 0.0)
                    {
                        infeasConstraint.Add(i);
                    }
                    else
                    {
                        normVal = 1.0;
                    }
                }

                if (normVal <= TOL_ABS)
                {
                    normVal = 1.0;
                }

                normA[i] = normVal;
                for (int j = 0; j < nvar; j++)
                {
                    A[i][j] /= normVal;
                }
                c[i] /= normVal;
            }

            if (infeasConstraint.Count > 0)
            {
                throw new Exception($"quadprog: unsolvable constraint with zero gradient. Indices: {string.Join(",", infeasConstraint)}");
            }

            tableau = Array.Empty<double[]>();
            if (nequ > 0)
            {
                double[][] Aequ = CreateMatrix(nequ, nvar, 0.0);
                double[] bequ = new double[nequ];
                for (int i = 0; i < nequ; i++)
                {
                    Aequ[i] = (double[])A[i].Clone();
                    bequ[i] = c[i];
                }

                int sort = 1;
                int act_cnt = QRfactorize(Aequ, out double[][] Qeq, out double[][] Req, out int[] perm, sort, tightToleranceRel);

                int rankQR = act_cnt;
                if (nequ > rankQR)
                {
                    double[][] Qcopy = CreateMatrix(Qeq.Length, Qeq[0].Length, 0.0);
                    for (int i = 0; i < Qeq.Length; i++)
                    {
                        Array.Copy(Qeq[i], Qcopy[i], Qeq[i].Length);
                    }
                    double[][] Rcopy = CreateMatrix(nequ, rankQR, 0.0);
                    for (int i = 0; i < nequ; i++)
                    {
                        for (int j = 0; j < rankQR; j++)
                        {
                            Rcopy[i][j] = Req[i][j];
                        }
                    }

                    bool append = QRappend(ref Qcopy, ref Rcopy, bequ, tightToleranceRel) != 0;
                    if (append)
                    {
                        var linDep = new bool[nequ];
                        for (int i = rankQR; i < nequ; i++)
                        {
                            for (int j = 0; j < nequ; j++)
                            {
                                if (Qeq[i][j] > tightToleranceAbs)
                                {
                                    linDep[j] = true;
                                }
                            }
                        }

                        infeasConstraint.Clear();
                        for (int i = 0; i < nequ; i++)
                        {
                            if (linDep[i])
                            {
                                infeasConstraint.Add(i);
                            }
                        }

                        throw new Exception($"quadprog: the equality constraints have no solution. Indices: {string.Join(",", infeasConstraint)}");
                    }
                }

                if (rankQR == nvar)
                {
                    double[] rhs = vecdotmat(Qeq, bequ);
                    QRbacksubstitute(out x, Req, perm, rankQR, rhs);

                    double[] tmp = matdotvecplusvec(A, x, c, -1.0);
                    infeasConstraint.Clear();
                    for (int i = nequ; i < ncst; i++)
                    {
                        if (tmp[i] > TOL_ABS)
                        {
                            infeasConstraint.Add(i);
                        }
                    }

                    if (infeasConstraint.Count > 0)
                    {
                        throw new Exception("quadprog: there is no feasible point; the equality constraints have been met");
                    }

                    actSet = new int[neqold + 2 * noFixed];
                    for (int i = 0; i < neqold; i++)
                    {
                        actSet[i] = i;
                    }
                    for (int i = 0; i < noFixed; i++)
                    {
                        actSet[neqold + i] = nold + fixInd[i];
                        actSet[neqold + noFixed + i] = nold + nvar + fixInd[i];
                    }

                    lambda = new double[nequ + nieq + lbdInd.Count + ubdInd.Count];
                    return -1;
                }
            }

            double[] residual = matdotvecplusvec(A, x, c, -1.0);
            int infeasCount = 0;
            for (int i = 0; i < nequ; i++)
            {
                if (Math.Abs(residual[i]) > TOL_ABS)
                {
                    infeasCount++;
                }
            }
            for (int i = nequ; i < ncst; i++)
            {
                if (residual[i] > TOL_ABS)
                {
                    infeasCount++;
                }
            }

            if (infeasCount > 0)
            {
                if (verbosity > 0)
                {
                    Console.WriteLine("starting point not feasible; use simplex");
                }

                int[] iposv;
                int[] izrov;
                tableau = CreateMatrix(ncst + 2, 2 * nvar + 2, 0.0);
                tableau[0][2 * nvar + 1] = -1.0;

                int info = 0;
                int k = 1;
                int mn = 2 * nvar + 2;
                int[] iposTemp = new int[ncst];
                int negCount = 0;
                for (int i = nequ; i < ncst; i++)
                {
                    if (c[i] >= 0.0)
                    {
                        tableau[k][mn - 1] = 1.0;
                        tableau[k][0] = c[i];
                        for (int j = 0; j < nvar; j++)
                        {
                            double value = A[i][j];
                            tableau[k][j + 1] = -value;
                            tableau[k][j + nvar + 1] = value;
                        }
                        k++;
                    }
                    else
                    {
                        iposTemp[negCount++] = i + 1;
                    }
                }

                mn = ncst - nequ - negCount + 1;
                for (int i = 0; i < negCount; i++)
                {
                    int pos = iposTemp[i] - 1;
                    tableau[mn + i][2 * nvar + 1] = -1.0;
                    tableau[mn + i][0] = -c[pos];
                    for (int j = 0; j < nvar; j++)
                    {
                        tableau[mn + i][j + 1] = A[pos][j];
                        tableau[mn + i][j + nvar + 1] = -A[pos][j];
                    }
                }

                mn = ncst - nequ + 1;
                for (int i = 0; i < nequ; i++)
                {
                    double sign = c[i] >= 0.0 ? 1.0 : -1.0;
                    tableau[i + mn][0] = sign * c[i];
                    for (int j = 0; j < nvar; j++)
                    {
                        tableau[i + mn][j + 1] = -sign * A[i][j];
                        tableau[i + mn][j + nvar + 1] = sign * A[i][j];
                    }
                }

                simplex(tableau, ncst - nequ - negCount, negCount, nequ, ref info, out izrov, out iposv, TOL_ABS);

                double gamma = 0.0;
                x = new double[nvar];
                for (int j = 1; j <= ncst; j++)
                {
                    if (iposv[j - 1] <= nvar)
                    {
                        x[iposv[j - 1] - 1] += tableau[j][0];
                    }
                    else if (iposv[j - 1] <= 2 * nvar)
                    {
                        x[iposv[j - 1] - nvar - 1] -= tableau[j][0];
                    }
                    else if (iposv[j - 1] == 2 * nvar + 1)
                    {
                        gamma = tableau[j][0];
                    }
                }

                if (info != 0 || gamma >= TOL_ABS)
                {
                    residual = matdotvecplusvec(A, x, c, -1.0);
                    infeasConstraint.Clear();
                    for (int i = 0; i < nequ; i++)
                    {
                        if (Math.Abs(residual[i]) > TOL_ABS)
                        {
                            infeasConstraint.Add(i);
                        }
                    }
                    for (int i = nequ; i < ncst; i++)
                    {
                        if (residual[i] > TOL_ABS)
                        {
                            infeasConstraint.Add(i);
                        }
                    }

                    if (infeasConstraint.Count > 0)
                    {
                        throw new Exception("quadprog(phase I): simplex optimum is not feasible");
                    }
                }
            }

            double[] gf = matdotvecplusvec(H, x, f);

            double[] alambda = Array.Empty<double>();
            double[] tlambda = Array.Empty<double>();

            int actCnt = 0;
            int[] actInd = new int[ncst + nvar];
            int actPtr = nequ;

            double[][] actTrn = Array.Empty<double[]>();
            double[][] Q = Array.Empty<double[]>();
            double[][] R = Array.Empty<double[]>();
            int rankQr = 0;
            double[][] Z;
            int[] solv = Array.Empty<int>();

            if (nequ > 0)
            {
                actTrn = CreateMatrix(nvar, nequ, 0.0);
                for (int i = 0; i < nvar; i++)
                {
                    for (int j = 0; j < nequ; j++)
                    {
                        actTrn[i][j] = A[j][i];
                    }
                }

                int sort = 1;
                actCnt = rankQr = QRfactorize(actTrn, out Q, out R, out solv, sort, tightToleranceRel);
                for (int i = 0; i < Q.Length; i++)
                {
                    for (int j = 0; j < Q[i].Length; j++)
                    {
                        if (Math.Abs(Q[i][j]) < MINDOUBLE)
                        {
                            Q[i][j] = 0.0;
                        }
                    }
                }
                for (int i = 0; i < R.Length; i++)
                {
                    for (int j = 0; j < R[i].Length; j++)
                    {
                        if (Math.Abs(R[i][j]) < MINDOUBLE)
                        {
                            R[i][j] = 0.0;
                        }
                    }
                }

                bool permutated = false;
                for (int i = 0; i < actCnt; i++)
                {
                    int pi = solv[i];
                    actInd[i] = pi;
                    if (pi != i && !permutated)
                    {
                        permutated = true;
                    }
                }

                if (actCnt < nequ || permutated)
                {
                    actTrn = CreateMatrix(nvar, actCnt, 0.0);
                    for (int i = 0; i < nvar; i++)
                    {
                        for (int j = 0; j < actCnt; j++)
                        {
                            actTrn[i][j] = A[actInd[j]][i];
                        }
                    }
                }

                Z = CreateMatrix(nvar, nvar, 0.0);
                for (int i = 0; i < nvar; i++)
                {
                    for (int j = actCnt; j < nvar; j++)
                    {
                        Z[i][j - actCnt] = Q[i][j];
                    }
                    for (int j = 0; j < actCnt; j++)
                    {
                        Z[i][nvar - actCnt + j] = 0.0;
                    }
                }
            }
            else
            {
                Z = CreateMatrix(nvar, nvar, 0.0);
                Q = CreateMatrix(nvar, nvar, 0.0);
                for (int i = 0; i < nvar; i++)
                {
                    for (int j = 0; j < nvar; j++)
                    {
                        double value = i == j ? 1.0 : 0.0;
                        Z[i][j] = value;
                        Q[i][j] = value;
                    }
                }
                solv = Array.Empty<int>();
                R = CreateMatrix(nvar, 0, 0.0);
            }

            double[][] Gz = projectMatrix(H, Z);
            double[] gz = vecdotmat(Z, gf);
            double[] tmpDir = Array.Empty<double>();
            int[] permTmp;
            QRsolve(ref tmpDir, out permTmp, ref Gz, gz, tightToleranceRel);
            double[] SD = matdotvec(Z, tmpDir);
            for (int i = 0; i < nvar; i++)
            {
                SD[i] = -SD[i];
            }

            double[] cstr = matdotvecplusvec(A, x, c, -1.0);
            int olin = 0;

            var gsd = new double[ncst];
            var gsdInd = new List<int>();
            var gsdIn2 = new List<int>();
            var dist = new double[0];
            var lamInd = new int[ncst];
            var lamInd2 = new int[ncst];

            for (int iters = 0; iters < MAXITER; iters++)
            {
                double[] gsdVec = matdotvec(A, SD);
                Array.Copy(gsdVec, gsd, ncst);
                double gsdNorm = norm(gsd);
                double tol = TOL_REL * gsdNorm;

                gsdInd.Clear();
                for (int i = 0; i < ncst; i++)
                {
                    bool flag = false;
                    for (int j = 0; j < actCnt; j++)
                    {
                        if (actInd[j] == i)
                        {
                            flag = true;
                            break;
                        }
                    }
                    if (!flag && gsd[i] > tol)
                    {
                        gsdInd.Add(i);
                    }
                }

                gsdIn2.Clear();
                double dmin = BIG;
                if (gsdInd.Count > 0)
                {
                    dist = new double[gsdInd.Count];
                    for (int i = 0; i < gsdInd.Count; i++)
                    {
                        int idx = gsdInd[i];
                        double dst = Math.Max(0.0, -cstr[idx] / gsd[idx]);
                        dist[i] = dst;
                        if (dst < dmin)
                        {
                            dmin = dst;
                        }
                    }
                }

                if (dmin < 1.0)
                {
                    for (int i = 0; i < gsdInd.Count; i++)
                    {
                        if (Math.Abs(dist[i] - dmin) <= TOL_REL * dmin)
                        {
                            gsdIn2.Add(gsdInd[i]);
                        }
                    }
                }

                int ind = -1;
                if (gsdIn2.Count > 0)
                {
                    ind = gsdIn2[0];
                    for (int j = 1; j < gsdIn2.Count; j++)
                    {
                        if (gsdIn2[j] < ind)
                        {
                            ind = gsdIn2[j];
                        }
                    }
                }

                if (dmin >= 1.0)
                {
                    for (int i = 0; i < nvar; i++)
                    {
                        x[i] += SD[i];
                    }

                    if (actCnt > 0)
                    {
                        double[] tmp = matdotvecplusvec(H, x, f);
                        double[] tmp2 = vecdotmat(Q, tmp);
                        for (int i = 0; i < nvar; i++)
                        {
                            tmp2[i] = -tmp2[i];
                        }

                        tlambda = new double[actCnt];
                        QRbacksubstitute(out tlambda, R, solv, rankQr, tmp2);

                        alambda = new double[actCnt];
                        for (int i = 0; i < nequ; i++)
                        {
                            alambda[i] = Math.Abs(tlambda[i]);
                        }
                        for (int i = nequ; i < actCnt; i++)
                        {
                            alambda[i] = tlambda[i];
                        }

                        int lamCnt = 0;
                        for (int i = 0; i < actCnt; i++)
                        {
                            if (alambda[i] < -TOL_ABS)
                            {
                                lamInd[lamCnt] = actInd[i];
                                lamInd2[lamCnt] = i;
                                lamCnt++;
                            }
                        }

                        if (lamCnt == 0)
                        {
                            actSet = new int[actCnt + noFixed];
                            lambda = new double[nequ + nieq + lbdInd.Count + ubdInd.Count];
                            for (int i = 0; i < actCnt; i++)
                            {
                                int ii = actInd[i];
                                double lmbd = tlambda[i] / normA[ii];
                                if (ii < neqold)
                                {
                                    lambda[ii] = lmbd;
                                    actSet[i] = ii;
                                }
                                else if (ii < nequ)
                                {
                                    int xtrEq = ii - neqold;
                                    int fixedVar = fixInd[xtrEq];
                                    int trInd = neqold + nieq + fixedVar;
                                    actSet[i] = trInd;
                                    actSet[actCnt + xtrEq] = trInd + nvar;
                                    lambda[trInd] = lmbd;
                                    lambda[trInd + nvar] = lmbd;
                                }
                                else if (ii < nequ + nieq)
                                {
                                    int trInd = ii - noFixed;
                                    actSet[i] = trInd;
                                    lambda[trInd] = lmbd;
                                }
                                else if (ii < nequ + nieq + lbdInd.Count)
                                {
                                    int idxLbd = ii - nequ - nieq;
                                    int trInd = neqold + nieq + lbdInd[idxLbd];
                                    actSet[i] = trInd;
                                    lambda[trInd] = lmbd;
                                }
                                else
                                {
                                    int idxUbd = ii - nequ - nieq - lbdInd.Count;
                                    int trInd = neqold + nieq + nvar + ubdInd[idxUbd];
                                    actSet[i] = trInd;
                                    lambda[trInd] = lmbd;
                                }
                            }

                            return iters;
                        }

                        int lind = ncst;
                        int k = 0;
                        for (int i = 0; i < lamCnt; i++)
                        {
                            if (lamInd[i] < lind)
                            {
                                lind = lamInd[i];
                                k = lamInd2[i];
                            }
                        }

                        for (int i = k; i < actCnt - 1; i++)
                        {
                            actInd[i] = actInd[i + 1];
                        }
                        actCnt--;

                        if (actCnt > 0)
                        {
                            actTrn = CreateMatrix(nvar, actCnt, 0.0);
                            for (int i = 0; i < nvar; i++)
                            {
                                for (int j = 0; j < actCnt; j++)
                                {
                                    actTrn[i][j] = A[actInd[j]][i];
                                }
                            }
                        }

                        int Acol = solv[k];
                        QRdelete(ref Q, ref R, Acol, tightToleranceAbs);
                        rankQr--;

                        for (int i = k; i < actCnt; i++)
                        {
                            solv[i] = solv[i + 1];
                        }
                        Array.Resize(ref solv, actCnt);
                        for (int i = 0; i < actCnt; i++)
                        {
                            if (solv[i] > Acol)
                            {
                                solv[i] -= 1;
                            }
                        }
                    }
                    else
                    {
                        actSet = new int[actCnt];
                        lambda = new double[nequ + nieq + lbdInd.Count + ubdInd.Count];
                        return iters;
                    }
                }
                else
                {
                    for (int i = 0; i < nvar; i++)
                    {
                        SD[i] *= dmin;
                        x[i] += SD[i];
                    }
                }

                gf = matdotvecplusvec(H, x, f);
                double[] tmpUpdate = matdotvec(A, SD);
                for (int i = 0; i < ncst; i++)
                {
                    cstr[i] += tmpUpdate[i];
                }

                int violated = 0;
                double tol2 = TOL_ABS;
                for (int i = 0; i < nequ; i++)
                {
                    if (Math.Abs(cstr[i]) > tol2)
                    {
                        violated++;
                        break;
                    }
                }
                for (int i = nequ; i < ncst; i++)
                {
                    if (cstr[i] > tol2)
                    {
                        violated++;
                        break;
                    }
                }

                if (violated > 0)
                {
                    throw new Exception("quadprog: the problem is badly conditioned");
                }

                if (ind >= 0 && dmin < 1.0)
                {
                    actPtr = actCnt;
                    actInd[actPtr] = ind;
                    actCnt++;

                    if (actTrn.Length == 0)
                    {
                        actTrn = CreateMatrix(nvar, actCnt, 0.0);
                    }
                    else
                    {
                        actTrn = CreateMatrix(nvar, actCnt, 0.0);
                    }

                    double[] newCnt = new double[nvar];
                    for (int i = 0; i < nvar; i++)
                    {
                        newCnt[i] = A[ind][i];
                        actTrn[i][actPtr] = newCnt[i];
                    }

                    if (Q.Length == 0)
                    {
                        Q = CreateMatrix(nvar, nvar, 0.0);
                        for (int i = 0; i < nvar; i++)
                        {
                            Q[i][i] = 1.0;
                        }
                    }

                    if (R.Length == 0)
                    {
                        R = CreateMatrix(nvar, 0, 0.0);
                    }

                    if (QRappend(ref Q, ref R, newCnt, tightToleranceRel) != 0)
                    {
                        Array.Resize(ref solv, actCnt);
                        solv[actCnt - 1] = actCnt - 1;
                        rankQr++;
                    }
                }

                if (actCnt < nvar)
                {
                    Z = CreateMatrix(nvar, nvar - actCnt, 0.0);
                    for (int i = 0; i < nvar; i++)
                    {
                        for (int j = actCnt; j < nvar; j++)
                        {
                            Z[i][j - actCnt] = Q[i][j];
                        }
                    }
                    actPtr = actCnt;
                    olin = 0;
                }
                else
                {
                    double[] tmpLambda = vecdotmat(Q, gf);
                    for (int i = 0; i < nvar; i++)
                    {
                        tmpLambda[i] = -tmpLambda[i];
                    }

                    tlambda = new double[actCnt];
                    QRbacksubstitute(out tlambda, R, solv, rankQr, tmpLambda);

                    alambda = new double[actCnt];
                    for (int i = 0; i < nequ; i++)
                    {
                        alambda[i] = Math.Abs(alambda[i]);
                    }
                    for (int i = nequ; i < actCnt; i++)
                    {
                        alambda[i] = tlambda[i];
                    }

                    int lamCnt = 0;
                    for (int i = 0; i < actCnt; i++)
                    {
                        if (alambda[i] < -TOL_ABS)
                        {
                            lamInd[lamCnt] = actInd[i];
                            lamInd2[lamCnt] = i;
                            lamCnt++;
                        }
                    }

                    if (lamCnt == 0)
                    {
                        actSet = new int[actCnt + noFixed];
                        lambda = new double[nequ + nieq + lbdInd.Count + ubdInd.Count];
                        for (int i = 0; i < actCnt; i++)
                        {
                            int ii = actInd[i];
                            double lmbd = tlambda[i] / normA[ii];
                            if (ii < neqold)
                            {
                                lambda[ii] = lmbd;
                                actSet[i] = ii;
                            }
                            else if (ii < nequ)
                            {
                                int xtrEq = ii - neqold;
                                int fixedVar = fixInd[xtrEq];
                                int trInd = neqold + nieq + fixedVar;
                                actSet[i] = trInd;
                                actSet[actCnt + xtrEq] = trInd + nvar;
                                lambda[trInd] = lmbd;
                                lambda[trInd + nvar] = lmbd;
                            }
                            else if (ii < nequ + nieq)
                            {
                                int trInd = ii - noFixed;
                                actSet[i] = trInd;
                                lambda[trInd] = lmbd;
                            }
                            else if (ii < nequ + nieq + lbdInd.Count)
                            {
                                int idxLbd = ii - nequ - nieq;
                                int trInd = neqold + nieq + lbdInd[idxLbd];
                                actSet[i] = trInd;
                                lambda[trInd] = lmbd;
                            }
                            else
                            {
                                int idxUbd = ii - nequ - nieq - lbdInd.Count;
                                int trInd = neqold + nieq + nvar + ubdInd[idxUbd];
                                actSet[i] = trInd;
                                lambda[trInd] = lmbd;
                            }
                        }

                        return iters;
                    }

                    int lind = ncst;
                    int k = 0;
                    for (int i = 0; i < lamCnt; i++)
                    {
                        if (lamInd[i] < lind)
                        {
                            lind = lamInd[i];
                            k = lamInd2[i];
                        }
                    }

                    for (int i = k; i < actCnt - 1; i++)
                    {
                        actInd[i] = actInd[i + 1];
                    }
                    actCnt--;

                    if (actCnt > 0)
                    {
                        actTrn = CreateMatrix(nvar, actCnt, 0.0);
                        for (int i = 0; i < nvar; i++)
                        {
                            for (int j = 0; j < actCnt; j++)
                            {
                                actTrn[i][j] = A[actInd[j]][i];
                            }
                        }
                    }

                    int Acol = solv[k];
                    QRdelete(ref Q, ref R, Acol, tightToleranceAbs);
                    rankQr--;

                    for (int i = k; i < actCnt; i++)
                    {
                        solv[i] = solv[i + 1];
                    }
                    Array.Resize(ref solv, actCnt);
                    for (int i = 0; i < actCnt; i++)
                    {
                        if (solv[i] > Acol)
                        {
                            solv[i] -= 1;
                        }
                    }

                    Z = CreateMatrix(nvar, nvar - actCnt, 0.0);
                    for (int i = 0; i < nvar; i++)
                    {
                        for (int j = actCnt; j < nvar; j++)
                        {
                            Z[i][j - actCnt] = Q[i][j];
                        }
                    }

                    olin = actInd[actPtr] + 1;
                }

                if (actCnt < nvar)
                {
                    gz = vecdotmat(Z, gf);
                    double normGz = norm(gz);
                    if (normGz < TOL_ABS)
                    {
                        SD = new double[nvar];
                    }
                    else
                    {
                        Gz = projectMatrix(H, Z);
                        int[] perm;
                        double[] diagR;
                        int defect2 = QRsolveWithDiagR(out tmpDir, out perm, out diagR, ref Gz, gz, tightToleranceRel);
                        SD = matdotvec(Z, tmpDir);
                        for (int i = 0; i < nvar; i++)
                        {
                            SD[i] = -SD[i];
                        }
                        if (verbosity > 1)
                        {
                            Console.WriteLine($"search direction: defect(Gz)= {defect2}, rank(Gz)= {diagR.Length}");
                        }
                    }
                }

                if (actCnt == nvar && olin != 0)
                {
                    for (int i = actPtr; i < actCnt - 1; i++)
                    {
                        actInd[i] = actInd[i + 1];
                    }
                    actPtr = nvar - 1;
                }
            }

            throw new Exception("quadprog: maximum number of iterations exceeded");
        }

        private static double[][] CreateMatrix(int rows, int cols, double value)
        {
            return new double[0][].Resize(rows, cols, value);
        }
    }
}
