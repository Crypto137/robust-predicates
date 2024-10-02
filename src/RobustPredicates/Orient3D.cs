using System;

namespace RobustPredicates
{
    public class Orient3D
    {
        public static unsafe double Fast(double* pa, double* pb, double* pc, double* pd)
        {
            double adx = pa[0] - pd[0];
            double bdx = pb[0] - pd[0];
            double cdx = pc[0] - pd[0];
            double ady = pa[1] - pd[1];
            double bdy = pb[1] - pd[1];
            double cdy = pc[1] - pd[1];
            double adz = pa[2] - pd[2];
            double bdz = pb[2] - pd[2];
            double cdz = pc[2] - pd[2];

            return adx * (bdy * cdz - bdz * cdy)
                 + bdx * (cdy * adz - cdz * ady)
                 + cdx * (ady * bdz - adz * bdy);
        }

        public static unsafe double Fast(ReadOnlySpan<double> pa, ReadOnlySpan<double> pb, ReadOnlySpan<double> pc, ReadOnlySpan<double> pd)
        {
            fixed (double* paPtr = pa)
            fixed (double* pbPtr = pb)
            fixed (double* pcPtr = pc)
            fixed (double* pdPtr = pd)
            {
                return Fast(paPtr, pbPtr, pcPtr, pdPtr);
            }
        }

        public static unsafe double Robust(double* pa, double* pb, double* pc, double* pd)
        {
            double adx = pa[0] - pd[0];
            double bdx = pb[0] - pd[0];
            double cdx = pc[0] - pd[0];
            double ady = pa[1] - pd[1];
            double bdy = pb[1] - pd[1];
            double cdy = pc[1] - pd[1];
            double adz = pa[2] - pd[2];
            double bdz = pb[2] - pd[2];
            double cdz = pc[2] - pd[2];

            double bdxcdy = bdx * cdy;
            double cdxbdy = cdx * bdy;

            double cdxady = cdx * ady;
            double adxcdy = adx * cdy;

            double adxbdy = adx * bdy;
            double bdxady = bdx * ady;

            double det = adz * (bdxcdy - cdxbdy)
                 + bdz * (cdxady - adxcdy)
                 + cdz * (adxbdy - bdxady);

            double permanent = (Math.Abs(bdxcdy) + Math.Abs(cdxbdy)) * Math.Abs(adz)
                       + (Math.Abs(cdxady) + Math.Abs(adxcdy)) * Math.Abs(bdz)
                       + (Math.Abs(adxbdy) + Math.Abs(bdxady)) * Math.Abs(cdz);

            double errbound = MacrosHelpers.O3derrboundA * permanent;
            if ((det > errbound) || (-det > errbound))
            {
                return det;
            }

            return Adapt(pa, pb, pc, pd, permanent);
        }

        public static unsafe double Robust(ReadOnlySpan<double> pa, ReadOnlySpan<double> pb, ReadOnlySpan<double> pc, ReadOnlySpan<double> pd)
        {
            fixed (double* paPtr = pa)
            fixed (double* pbPtr = pb)
            fixed (double* pcPtr = pc)
            fixed (double* pdPtr = pd)
            {
                return Robust(paPtr, pbPtr, pcPtr, pdPtr);
            }
        }

        private static unsafe double Adapt(double* pa, double* pb, double* pc, double* pd, double permanent)
        {
            double adx = pa[0] - pd[0];
            double bdx = pb[0] - pd[0];
            double cdx = pc[0] - pd[0];
            double ady = pa[1] - pd[1];
            double bdy = pb[1] - pd[1];
            double cdy = pc[1] - pd[1];
            double adz = pa[2] - pd[2];
            double bdz = pb[2] - pd[2];
            double cdz = pc[2] - pd[2];

            MacrosHelpers.TwoProduct(bdx, cdy, out double bdxcdy1, out double bdxcdy0);
            MacrosHelpers.TwoProduct(cdx, bdy, out double cdxbdy1, out double cdxbdy0);
            double* bc = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, out bc[3], out bc[2], out bc[1], out bc[0]);
            double* adet = stackalloc double[8];
            int alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bc, adz, adet);

            MacrosHelpers.TwoProduct(cdx, ady, out double cdxady1, out double cdxady0);
            MacrosHelpers.TwoProduct(adx, cdy, out double adxcdy1, out double adxcdy0);
            double* ca = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(cdxady1, cdxady0, adxcdy1, adxcdy0, out ca[3], out ca[2], out ca[1], out ca[0]);
            double* bdet = stackalloc double[8];
            int blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ca, bdz, bdet);

            MacrosHelpers.TwoProduct(adx, bdy, out double adxbdy1, out double adxbdy0);
            MacrosHelpers.TwoProduct(bdx, ady, out double bdxady1, out double bdxady0);
            double* ab = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(adxbdy1, adxbdy0, bdxady1, bdxady0, out ab[3], out ab[2], out ab[1], out ab[0]);

            double* cdet = stackalloc double[8];
            int clen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ab, cdz, cdet);

            double* abdet = stackalloc double[16];
            int ablen = ArithmeticFunctions.FastExpansionSumZeroelim(alen, adet, blen, bdet, abdet);
            double* fin1 = stackalloc double[192];
            int finlength = ArithmeticFunctions.FastExpansionSumZeroelim(ablen, abdet, clen, cdet, fin1);

            double det = ArithmeticFunctions.Estimate(finlength, fin1);

            double errbound = MacrosHelpers.O3derrboundB * permanent;
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            MacrosHelpers.TwoDiffTail(pa[0], pd[0], adx, out double adxtail);
            MacrosHelpers.TwoDiffTail(pb[0], pd[0], bdx, out double bdxtail);
            MacrosHelpers.TwoDiffTail(pc[0], pd[0], cdx, out double cdxtail);
            MacrosHelpers.TwoDiffTail(pa[1], pd[1], ady, out double adytail);
            MacrosHelpers.TwoDiffTail(pb[1], pd[1], bdy, out double bdytail);
            MacrosHelpers.TwoDiffTail(pc[1], pd[1], cdy, out double cdytail);
            MacrosHelpers.TwoDiffTail(pa[2], pd[2], adz, out double adztail);
            MacrosHelpers.TwoDiffTail(pb[2], pd[2], bdz, out double bdztail);
            MacrosHelpers.TwoDiffTail(pc[2], pd[2], cdz, out double cdztail);

            if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
                && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)
                && (adztail == 0.0) && (bdztail == 0.0) && (cdztail == 0.0))
            {
                return det;
            }

            errbound = MacrosHelpers.O3derrboundC * permanent + MacrosHelpers.Resulterrbound * Math.Abs(det);
            det += (adz * ((bdx * cdytail + cdy * bdxtail)
                           - (bdy * cdxtail + cdx * bdytail))
                    + adztail * (bdx * cdy - bdy * cdx))
                 + (bdz * ((cdx * adytail + ady * cdxtail)
                           - (cdy * adxtail + adx * cdytail))
                    + bdztail * (cdx * ady - cdy * adx))
                 + (cdz * ((adx * bdytail + bdy * adxtail)
                           - (ady * bdxtail + bdx * adytail))
                    + cdztail * (adx * bdy - ady * bdx));
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            double* finnow = fin1;
            double* fin2 = stackalloc double[192];
            double* finother = fin2;

            double* at_b = stackalloc double[4];
            double* at_c = stackalloc double[4];
            int at_blen;
            int at_clen;
            if (adxtail == 0.0)
            {
                if (adytail == 0.0)
                {
                    at_b[0] = 0.0;
                    at_blen = 1;
                    at_c[0] = 0.0;
                    at_clen = 1;
                }
                else
                {
                    double negate = -adytail;
                    MacrosHelpers.TwoProduct(negate, bdx, out at_b[1], out at_b[0]);
                    at_blen = 2;
                    MacrosHelpers.TwoProduct(adytail, cdx, out at_c[1], out at_c[0]);
                    at_clen = 2;
                }
            }
            else
            {
                if (adytail == 0.0)
                {
                    MacrosHelpers.TwoProduct(adxtail, bdy, out at_b[1], out at_b[0]);
                    at_blen = 2;
                    double negate = -adxtail;
                    MacrosHelpers.TwoProduct(negate, cdy, out at_c[1], out at_c[0]);
                    at_clen = 2;
                }
                else
                {
                    MacrosHelpers.TwoProduct(adxtail, bdy, out double adxt_bdy1, out double adxt_bdy0);
                    MacrosHelpers.TwoProduct(adytail, bdx, out double adyt_bdx1, out double adyt_bdx0);
                    MacrosHelpers.TwoTwoDiff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0,
                   out double at_blarge, out at_b[2], out at_b[1], out at_b[0]);
                    at_b[3] = at_blarge;
                    at_blen = 4;
                    MacrosHelpers.TwoProduct(adytail, cdx, out double adyt_cdx1, out double adyt_cdx0);
                    MacrosHelpers.TwoProduct(adxtail, cdy, out double adxt_cdy1, out double adxt_cdy0);
                    MacrosHelpers.TwoTwoDiff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0,
                   out at_c[3], out at_c[2], out at_c[1], out at_c[0]);
                    at_clen = 4;
                }
            }
            double* bt_c = stackalloc double[4];
            double* bt_a = stackalloc double[4];
            int bt_alen;
            int bt_clen;
            if (bdxtail == 0.0)
            {
                if (bdytail == 0.0)
                {
                    bt_c[0] = 0.0;
                    bt_clen = 1;
                    bt_a[0] = 0.0;
                    bt_alen = 1;
                }
                else
                {
                    double negate = -bdytail;
                    MacrosHelpers.TwoProduct(negate, cdx, out bt_c[1], out bt_c[0]);
                    bt_clen = 2;
                    MacrosHelpers.TwoProduct(bdytail, adx, out bt_a[1], out bt_a[0]);
                    bt_alen = 2;
                }
            }
            else
            {
                if (bdytail == 0.0)
                {
                    MacrosHelpers.TwoProduct(bdxtail, cdy, out bt_c[1], out bt_c[0]);
                    bt_clen = 2;
                    double negate = -bdxtail;
                    MacrosHelpers.TwoProduct(negate, ady, out bt_a[1], out bt_a[0]);
                    bt_alen = 2;
                }
                else
                {
                    MacrosHelpers.TwoProduct(bdxtail, cdy, out double bdxt_cdy1, out double bdxt_cdy0);
                    MacrosHelpers.TwoProduct(bdytail, cdx, out double bdyt_cdx1, out double bdyt_cdx0);
                    MacrosHelpers.TwoTwoDiff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0,
                    out bt_c[3], out bt_c[2], out bt_c[1], out bt_c[0]);
                    bt_clen = 4;
                    MacrosHelpers.TwoProduct(bdytail, adx, out double bdyt_adx1, out double bdyt_adx0);
                    MacrosHelpers.TwoProduct(bdxtail, ady, out double bdxt_ady1, out double bdxt_ady0);
                    MacrosHelpers.TwoTwoDiff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0,
                 out bt_a[3], out bt_a[2], out bt_a[1], out bt_a[0]);
                    bt_alen = 4;
                }
            }
            double* ct_a = stackalloc double[4];
            double* ct_b = stackalloc double[4];
            int ct_alen;
            int ct_blen;
            if (cdxtail == 0.0)
            {
                if (cdytail == 0.0)
                {
                    ct_a[0] = 0.0;
                    ct_alen = 1;
                    ct_b[0] = 0.0;
                    ct_blen = 1;
                }
                else
                {
                    double negate = -cdytail;
                    MacrosHelpers.TwoProduct(negate, adx, out ct_a[1], out ct_a[0]);
                    ct_alen = 2;
                    MacrosHelpers.TwoProduct(cdytail, bdx, out ct_b[1], out ct_b[0]);
                    ct_blen = 2;
                }
            }
            else
            {
                if (cdytail == 0.0)
                {
                    MacrosHelpers.TwoProduct(cdxtail, ady, out ct_a[1], out ct_a[0]);
                    ct_alen = 2;
                    double negate = -cdxtail;
                    MacrosHelpers.TwoProduct(negate, bdy, out ct_b[1], out ct_b[0]);
                    ct_blen = 2;
                }
                else
                {
                    MacrosHelpers.TwoProduct(cdxtail, ady, out double cdxt_ady1, out double cdxt_ady0);
                    MacrosHelpers.TwoProduct(cdytail, adx, out double cdyt_adx1, out double cdyt_adx0);
                    MacrosHelpers.TwoTwoDiff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0,
                   out ct_a[3], out ct_a[2], out ct_a[1], out ct_a[0]);
                    ct_alen = 4;
                    MacrosHelpers.TwoProduct(cdytail, bdx, out double cdyt_bdx1, out double cdyt_bdx0);
                    MacrosHelpers.TwoProduct(cdxtail, bdy, out double cdxt_bdy1, out double cdxt_bdy0);
                    MacrosHelpers.TwoTwoDiff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0,
                   out ct_b[3], out ct_b[2], out ct_b[1], out ct_b[0]);
                    ct_blen = 4;
                }
            }
            double* bct = stackalloc double[8];
            int bctlen = ArithmeticFunctions.FastExpansionSumZeroelim(bt_clen, bt_c, ct_blen, ct_b, bct);
            double* u = stackalloc double[4];
            double* v = stackalloc double[12];
            double* w = stackalloc double[16];
            int wlength = ArithmeticFunctions.ScaleExpansionZeroelim(bctlen, bct, adz, w);
            finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, wlength, w,
                                                    finother);
            double* finswap = finnow; finnow = finother; finother = finswap;
            double* cat = stackalloc double[8];
            int catlen = ArithmeticFunctions.FastExpansionSumZeroelim(ct_alen, ct_a, at_clen, at_c, cat);
            wlength = ArithmeticFunctions.ScaleExpansionZeroelim(catlen, cat, bdz, w);
            finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, wlength, w,
                                                    finother);
            finswap = finnow; finnow = finother; finother = finswap;
            double* abt = stackalloc double[8];
            int abtlen = ArithmeticFunctions.FastExpansionSumZeroelim(at_blen, at_b, bt_alen, bt_a, abt);
            wlength = ArithmeticFunctions.ScaleExpansionZeroelim(abtlen, abt, cdz, w);
            finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, wlength, w,
                                                    finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;

            if (adztail != 0.0)
            {
                int vlength = ArithmeticFunctions.ScaleExpansionZeroelim(4, bc, adztail, v);
                finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, vlength, v,
                                                        finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (bdztail != 0.0)
            {
                int vlength = ArithmeticFunctions.ScaleExpansionZeroelim(4, ca, bdztail, v);
                finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, vlength, v,
                                                        finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (cdztail != 0.0)
            {
                int vlength = ArithmeticFunctions.ScaleExpansionZeroelim(4, ab, cdztail, v);
                finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, vlength, v,
                                                        finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }

            if (adxtail != 0.0)
            {
                if (bdytail != 0.0)
                {
                    MacrosHelpers.TwoProduct(adxtail, bdytail, out double adxt_bdyt1, out double adxt_bdyt0);
                    MacrosHelpers.TwoOneProduct(adxt_bdyt1, adxt_bdyt0, cdz, out u[3], out u[2], out u[1], out u[0]);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, 4, u, finother);
                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;
                    if (cdztail != 0.0)
                    {
                        MacrosHelpers.TwoOneProduct(adxt_bdyt1, adxt_bdyt0, cdztail, out u[3], out u[2], out u[1], out u[0]);
                        finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, 4, u, finother);
                        finswap = finnow;
                        finnow = finother;
                        finother = finswap;
                    }
                }
                if (cdytail != 0.0)
                {
                    double negate = -adxtail;
                    MacrosHelpers.TwoProduct(negate, cdytail, out double adxt_cdyt1, out double adxt_cdyt0);
                    MacrosHelpers.TwoOneProduct(adxt_cdyt1, adxt_cdyt0, bdz, out u[3], out u[2], out u[1], out u[0]);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, 4, u, finother);
                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;
                    if (bdztail != 0.0)
                    {
                        MacrosHelpers.TwoOneProduct(adxt_cdyt1, adxt_cdyt0, bdztail, out u[3], out u[2], out u[1], out u[0]);
                        finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, 4, u, finother);
                        finswap = finnow; finnow = finother; finother = finswap;
                    }
                }
            }
            if (bdxtail != 0.0)
            {
                if (cdytail != 0.0)
                {
                    MacrosHelpers.TwoProduct(bdxtail, cdytail, out double bdxt_cdyt1, out double bdxt_cdyt0);
                    MacrosHelpers.TwoOneProduct(bdxt_cdyt1, bdxt_cdyt0, adz, out u[3], out u[2], out u[1], out u[0]);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, 4, u,
                                                            finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                    if (adztail != 0.0)
                    {
                        MacrosHelpers.TwoOneProduct(bdxt_cdyt1, bdxt_cdyt0, adztail, out u[3], out u[2], out u[1], out u[0]);
                        finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, 4, u,
                                                                finother);
                        finswap = finnow; finnow = finother; finother = finswap;
                    }
                }
                if (adytail != 0.0)
                {
                    double negate = -bdxtail;
                    MacrosHelpers.TwoProduct(negate, adytail, out double bdxt_adyt1, out double bdxt_adyt0);
                    MacrosHelpers.TwoOneProduct(bdxt_adyt1, bdxt_adyt0, cdz, out u[3], out u[2], out u[1], out u[0]);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, 4, u,
                                                            finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                    if (cdztail != 0.0)
                    {
                        MacrosHelpers.TwoOneProduct(bdxt_adyt1, bdxt_adyt0, cdztail, out u[3], out u[2], out u[1], out u[0]);
                        finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, 4, u,
                                                                finother);
                        finswap = finnow; finnow = finother; finother = finswap;
                    }
                }
            }
            if (cdxtail != 0.0)
            {
                if (adytail != 0.0)
                {
                    MacrosHelpers.TwoProduct(cdxtail, adytail, out double cdxt_adyt1, out double cdxt_adyt0);
                    MacrosHelpers.TwoOneProduct(cdxt_adyt1, cdxt_adyt0, bdz, out u[3], out u[2], out u[1], out u[0]);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, 4, u,
                                                            finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                    if (bdztail != 0.0)
                    {
                        MacrosHelpers.TwoOneProduct(cdxt_adyt1, cdxt_adyt0, bdztail, out u[3], out u[2], out u[1], out u[0]);
                        finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, 4, u,
                                                                finother);
                        finswap = finnow; finnow = finother; finother = finswap;
                    }
                }
                if (bdytail != 0.0)
                {
                    double negate = -cdxtail;
                    MacrosHelpers.TwoProduct(negate, bdytail, out double cdxt_bdyt1, out double cdxt_bdyt0);
                    MacrosHelpers.TwoOneProduct(cdxt_bdyt1, cdxt_bdyt0, adz, out u[3], out u[2], out u[1], out u[0]);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, 4, u,
                                                            finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                    if (adztail != 0.0)
                    {
                        MacrosHelpers.TwoOneProduct(cdxt_bdyt1, cdxt_bdyt0, adztail, out u[3], out u[2], out u[1], out u[0]);
                        finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, 4, u,
                                                                finother);
                        finswap = finnow; finnow = finother; finother = finswap;
                    }
                }
            }

            if (adztail != 0.0)
            {
                wlength = ArithmeticFunctions.ScaleExpansionZeroelim(bctlen, bct, adztail, w);
                finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, wlength, w,
                                                        finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (bdztail != 0.0)
            {
                wlength = ArithmeticFunctions.ScaleExpansionZeroelim(catlen, cat, bdztail, w);
                finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, wlength, w,
                                                        finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (cdztail != 0.0)
            {
                wlength = ArithmeticFunctions.ScaleExpansionZeroelim(abtlen, abt, cdztail, w);
                finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, wlength, w,
                                                        finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }

            return finnow[finlength - 1];
        }
    }
}
