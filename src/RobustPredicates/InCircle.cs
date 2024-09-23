using System;

namespace RobustPredicates
{
    public static class InCircle
    {
        public static double Fast(ReadOnlySpan<double> pa, ReadOnlySpan<double> pb, ReadOnlySpan<double> pc, ReadOnlySpan<double> pd)
        {
            double adx = pa[0] - pd[0];
            double ady = pa[1] - pd[1];
            double bdx = pb[0] - pd[0];
            double bdy = pb[1] - pd[1];
            double cdx = pc[0] - pd[0];
            double cdy = pc[1] - pd[1];

            double abdet = adx * bdy - bdx * ady;
            double bcdet = bdx * cdy - cdx * bdy;
            double cadet = cdx * ady - adx * cdy;
            double alift = adx * adx + ady * ady;
            double blift = bdx * bdx + bdy * bdy;
            double clift = cdx * cdx + cdy * cdy;

            return alift * bcdet + blift * cadet + clift * abdet;
        }

        internal static double Exact(ReadOnlySpan<double> pa, ReadOnlySpan<double> pb, ReadOnlySpan<double> pc, ReadOnlySpan<double> pd)
        {
            MacrosHelpers.TwoProduct(pa[0], pb[1], out double axby1, out double axby0);
            MacrosHelpers.TwoProduct(pb[0], pa[1], out double bxay1, out double bxay0);

            Span<double> ab = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(axby1, axby0, bxay1, bxay0, out ab[3], out ab[2], out ab[1], out ab[0]);

            Span<double> bc = stackalloc double[4];
            MacrosHelpers.TwoProduct(pb[0], pc[1], out double bxcy1, out double bxcy0);
            MacrosHelpers.TwoProduct(pc[0], pb[1], out double cxby1, out double cxby0);
            MacrosHelpers.TwoTwoDiff(bxcy1, bxcy0, cxby1, cxby0, out bc[3], out bc[2], out bc[1], out bc[0]);

            Span<double> cd = stackalloc double[4];
            MacrosHelpers.TwoProduct(pc[0], pd[1], out double cxdy1, out double cxdy0);
            MacrosHelpers.TwoProduct(pd[0], pc[1], out double dxcy1, out double dxcy0);
            MacrosHelpers.TwoTwoDiff(cxdy1, cxdy0, dxcy1, dxcy0, out cd[3], out cd[2], out cd[1], out cd[0]);

            Span<double> da = stackalloc double[4];
            MacrosHelpers.TwoProduct(pd[0], pa[1], out double dxay1, out double dxay0);
            MacrosHelpers.TwoProduct(pa[0], pd[1], out double axdy1, out double axdy0);
            MacrosHelpers.TwoTwoDiff(dxay1, dxay0, axdy1, axdy0, out da[3], out da[2], out da[1], out da[0]);

            Span<double> ac = stackalloc double[4];
            MacrosHelpers.TwoProduct(pa[0], pc[1], out double axcy1, out double axcy0);
            MacrosHelpers.TwoProduct(pc[0], pa[1], out double cxay1, out double cxay0);
            MacrosHelpers.TwoTwoDiff(axcy1, axcy0, cxay1, cxay0, out ac[3], out ac[2], out ac[1], out ac[0]);

            Span<double> bd = stackalloc double[4];
            MacrosHelpers.TwoProduct(pb[0], pd[1], out double bxdy1, out double bxdy0);
            MacrosHelpers.TwoProduct(pd[0], pb[1], out double dxby1, out double dxby0);
            MacrosHelpers.TwoTwoDiff(bxdy1, bxdy0, dxby1, dxby0, out bd[3], out bd[2], out bd[1], out bd[0]);

            Span<double> temp8 = stackalloc double[8];
            int templen = ArithmeticFunctions.FastExpansionSumZeroelim(4, cd, 4, da, temp8);
            Span<double> cda = stackalloc double[12];
            int cdalen = ArithmeticFunctions.FastExpansionSumZeroelim(templen, temp8, 4, ac, cda);
            Span<double> dab = stackalloc double[12];
            templen = ArithmeticFunctions.FastExpansionSumZeroelim(4, da, 4, ab, temp8);
            int dablen = ArithmeticFunctions.FastExpansionSumZeroelim(templen, temp8, 4, bd, dab);

            for (int i = 0; i < 4; i++)
            {
                bd[i] = -bd[i];
                ac[i] = -ac[i];
            }

            Span<double> abc = stackalloc double[12];
            templen = ArithmeticFunctions.FastExpansionSumZeroelim(4, ab, 4, bc, temp8);
            int abclen = ArithmeticFunctions.FastExpansionSumZeroelim(templen, temp8, 4, ac, abc);
            templen = ArithmeticFunctions.FastExpansionSumZeroelim(4, bc, 4, cd, temp8);

            Span<double> bcd = stackalloc double[12];
            int bcdlen = ArithmeticFunctions.FastExpansionSumZeroelim(templen, temp8, 4, bd, bcd);

            Span<double> det24x = stackalloc double[24];
            int xlen = ArithmeticFunctions.ScaleExpansionZeroelim(bcdlen, bcd, pa[0], det24x);
            Span<double> det48x = stackalloc double[48];
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(xlen, det24x, pa[0], det48x);
            Span<double> det24y = stackalloc double[24];
            int ylen = ArithmeticFunctions.ScaleExpansionZeroelim(bcdlen, bcd, pa[1], det24y);
            Span<double> det48y = stackalloc double[48];
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(ylen, det24y, pa[1], det48y);
            Span<double> adet = stackalloc double[96];
            int alen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, det48x, ylen, det48y, adet);

            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(cdalen, cda, pb[0], det24x);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(xlen, det24x, -pb[0], det48x);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(cdalen, cda, pb[1], det24y);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(ylen, det24y, -pb[1], det48y);
            Span<double> bdet = stackalloc double[96];
            int blen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, det48x, ylen, det48y, bdet);

            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(dablen, dab, pc[0], det24x);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(xlen, det24x, pc[0], det48x);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(dablen, dab, pc[1], det24y);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(ylen, det24y, pc[1], det48y);
            Span<double> cdet = stackalloc double[96];
            int clen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, det48x, ylen, det48y, cdet);

            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(abclen, abc, pd[0], det24x);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(xlen, det24x, -pd[0], det48x);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(abclen, abc, pd[1], det24y);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(ylen, det24y, -pd[1], det48y);
            Span<double> ddet = stackalloc double[96];
            int dlen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, det48x, ylen, det48y, ddet);

            Span<double> abdet = stackalloc double[192];
            int ablen = ArithmeticFunctions.FastExpansionSumZeroelim(alen, adet, blen, bdet, abdet);
            Span<double> cddet = stackalloc double[192];
            int cdlen = ArithmeticFunctions.FastExpansionSumZeroelim(clen, cdet, dlen, ddet, cddet);
            Span<double> deter = stackalloc double[384];
            int deterlen = ArithmeticFunctions.FastExpansionSumZeroelim(ablen, abdet, cdlen, cddet, deter);
            return deter[deterlen - 1];
        }

        private static double Adapt(ReadOnlySpan<double> pa, ReadOnlySpan<double> pb, ReadOnlySpan<double> pc, ReadOnlySpan<double> pd, double permanent)
        {
            double adx = pa[0] - pd[0];
            double bdx = pb[0] - pd[0];
            double cdx = pc[0] - pd[0];
            double ady = pa[1] - pd[1];
            double bdy = pb[1] - pd[1];
            double cdy = pc[1] - pd[1];

            MacrosHelpers.TwoProduct(bdx, cdy, out double bdxcdy1, out double bdxcdy0);
            MacrosHelpers.TwoProduct(cdx, bdy, out double cdxbdy1, out double cdxbdy0);
            Span<double> bc = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, out bc[3], out bc[2], out bc[1], out bc[0]);

            Span<double> axbc = stackalloc double[8];
            int axbclen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bc, adx, axbc);
            Span<double> axxbc = stackalloc double[16];
            int axxbclen = ArithmeticFunctions.ScaleExpansionZeroelim(axbclen, axbc, adx, axxbc);
            Span<double> aybc = stackalloc double[8];
            int aybclen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bc, ady, aybc);
            Span<double> ayybc = stackalloc double[16];
            int ayybclen = ArithmeticFunctions.ScaleExpansionZeroelim(aybclen, aybc, ady, ayybc);
            Span<double> adet = stackalloc double[32];
            int alen = ArithmeticFunctions.FastExpansionSumZeroelim(axxbclen, axxbc, ayybclen, ayybc, adet);

            MacrosHelpers.TwoProduct(cdx, ady, out double cdxady1, out double cdxady0);
            MacrosHelpers.TwoProduct(adx, cdy, out double adxcdy1, out double adxcdy0);
            Span<double> ca = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(cdxady1, cdxady0, adxcdy1, adxcdy0, out ca[3], out ca[2], out ca[1], out ca[0]);

            Span<double> bxca = stackalloc double[8];
            int bxcalen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ca, bdx, bxca);
            Span<double> bxxca = stackalloc double[16];
            int bxxcalen = ArithmeticFunctions.ScaleExpansionZeroelim(bxcalen, bxca, bdx, bxxca);
            Span<double> byca = stackalloc double[8];
            int bycalen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ca, bdy, byca);
            Span<double> byyca = stackalloc double[16];
            int byycalen = ArithmeticFunctions.ScaleExpansionZeroelim(bycalen, byca, bdy, byyca);
            Span<double> bdet = stackalloc double[32];
            int blen = ArithmeticFunctions.FastExpansionSumZeroelim(bxxcalen, bxxca, byycalen, byyca, bdet);



            MacrosHelpers.TwoProduct(adx, bdy, out double adxbdy1, out double adxbdy0);
            MacrosHelpers.TwoProduct(bdx, ady, out double bdxady1, out double bdxady0);
            Span<double> ab = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(adxbdy1, adxbdy0, bdxady1, bdxady0, out ab[3], out ab[2], out ab[1], out ab[0]);

            Span<double> cxab = stackalloc double[8];
            int cxablen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ab, cdx, cxab);
            Span<double> cxxab = stackalloc double[16];
            int cxxablen = ArithmeticFunctions.ScaleExpansionZeroelim(cxablen, cxab, cdx, cxxab);

            Span<double> cyab = stackalloc double[8];
            int cyablen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ab, cdy, cyab);
            Span<double> cyyab = stackalloc double[16];
            int cyyablen = ArithmeticFunctions.ScaleExpansionZeroelim(cyablen, cyab, cdy, cyyab);
            Span<double> cdet = stackalloc double[32];
            int clen = ArithmeticFunctions.FastExpansionSumZeroelim(cxxablen, cxxab, cyyablen, cyyab, cdet);

            Span<double> abdet = stackalloc double[64];
            int ablen = ArithmeticFunctions.FastExpansionSumZeroelim(alen, adet, blen, bdet, abdet);
            Span<double> fin1 = stackalloc double[1152];
            int finlength = ArithmeticFunctions.FastExpansionSumZeroelim(ablen, abdet, clen, cdet, fin1);

            double det = ArithmeticFunctions.Estimate(finlength, fin1);
            double errbound = MacrosHelpers.IccerrboundB * permanent;

            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            MacrosHelpers.TwoDiffTail(pa[0], pd[0], adx, out double adxtail);
            MacrosHelpers.TwoDiffTail(pa[1], pd[1], ady, out double adytail);
            MacrosHelpers.TwoDiffTail(pb[0], pd[0], bdx, out double bdxtail);
            MacrosHelpers.TwoDiffTail(pb[1], pd[1], bdy, out double bdytail);
            MacrosHelpers.TwoDiffTail(pc[0], pd[0], cdx, out double cdxtail);
            MacrosHelpers.TwoDiffTail(pc[1], pd[1], cdy, out double cdytail);

            if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
                && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0))
            {
                return det;
            }

            errbound = MacrosHelpers.IccerrboundC * permanent + MacrosHelpers.Resulterrbound * Math.Abs(det);
            det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail)
                                               - (bdy * cdxtail + cdx * bdytail))
                    + 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
                 + ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail)
                                               - (cdy * adxtail + adx * cdytail))
                    + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
                 + ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail)
                                               - (ady * bdxtail + bdx * adytail))
                    + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }


            Span<double> finnow = fin1;
            Span<double> fin2 = stackalloc double[1152];
            Span<double> finother = fin2;
            Span<double> aa = stackalloc double[4];
            Span<double> bb = stackalloc double[4];
            Span<double> cc = stackalloc double[4];

            if ((bdxtail != 0.0) || (bdytail != 0.0)
                || (cdxtail != 0.0) || (cdytail != 0.0))
            {
                MacrosHelpers.Square(adx, out double adxadx1, out double adxadx0);
                MacrosHelpers.Square(ady, out double adyady1, out double adyady0);
                MacrosHelpers.TwoTwoSum(adxadx1, adxadx0, adyady1, adyady0, out aa[3], out aa[2], out aa[1], out aa[0]);
            }
            if ((cdxtail != 0.0) || (cdytail != 0.0)
                || (adxtail != 0.0) || (adytail != 0.0))
            {
                MacrosHelpers.Square(bdx, out double bdxbdx1, out double bdxbdx0);
                MacrosHelpers.Square(bdy, out double bdybdy1, out double bdybdy0);
                MacrosHelpers.TwoTwoSum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, out bb[3], out bb[2], out bb[1], out bb[0]);
            }
            if ((adxtail != 0.0) || (adytail != 0.0)
                || (bdxtail != 0.0) || (bdytail != 0.0))
            {
                MacrosHelpers.Square(cdx, out double cdxcdx1, out double cdxcdx0);
                MacrosHelpers.Square(cdy, out double cdycdy1, out double cdycdy0);
                MacrosHelpers.TwoTwoSum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, out cc[3], out cc[2], out cc[1], out cc[0]);
            }


            Span<double> finswap = stackalloc double[0];
            Span<double> temp16a = stackalloc double[16];
            Span<double> temp32a = stackalloc double[32];
            Span<double> temp16b = stackalloc double[16];
            Span<double> temp32b = stackalloc double[32];
            Span<double> temp16c = stackalloc double[16];
            Span<double> temp48 = stackalloc double[48];
            Span<double> axtbc = stackalloc double[8];
            int axtbclen = 0;
            if (adxtail != 0.0)
            {
                axtbclen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bc, adxtail, axtbc);
                int temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(axtbclen, axtbc, 2.0 * adx, temp16a);
                Span<double> axtcc = stackalloc double[8];
                int axtcclen = ArithmeticFunctions.ScaleExpansionZeroelim(4, cc, adxtail, axtcc);
                int temp16blen = ArithmeticFunctions.ScaleExpansionZeroelim(axtcclen, axtcc, bdy, temp16b);
                Span<double> axtbb = stackalloc double[8];
                int axtbblen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bb, adxtail, axtbb);
                int temp16clen = ArithmeticFunctions.ScaleExpansionZeroelim(axtbblen, axtbb, -cdy, temp16c);

                int temp32alen = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                        temp16blen, temp16b, temp32a);

                int temp48len = ArithmeticFunctions.FastExpansionSumZeroelim(temp16clen, temp16c,
                                                        temp32alen, temp32a, temp48);
                finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                        temp48, finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap;
            }

            int aytbclen = 0;
            Span<double> aytbc = stackalloc double[8];
            if (adytail != 0.0)
            {
                aytbclen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bc, adytail, aytbc);
                int temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(aytbclen, aytbc, 2.0 * ady, temp16a);
                Span<double> aytbb = stackalloc double[8];
                int aytbblen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bb, adytail, aytbb);
                int temp16blen = ArithmeticFunctions.ScaleExpansionZeroelim(aytbblen, aytbb, cdx, temp16b);
                Span<double> aytcc = stackalloc double[8];
                int aytcclen = ArithmeticFunctions.ScaleExpansionZeroelim(4, cc, adytail, aytcc);
                int temp16clen = ArithmeticFunctions.ScaleExpansionZeroelim(aytcclen, aytcc, -bdx, temp16c);
                int temp32alen = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                        temp16blen, temp16b, temp32a);
                int temp48len = ArithmeticFunctions.FastExpansionSumZeroelim(temp16clen, temp16c,
                                                        temp32alen, temp32a, temp48);
                finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                        temp48, finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap;
            }

            int bxtcalen = 0;
            Span<double> bxtca = stackalloc double[8];
            if (bdxtail != 0.0)
            {
                bxtcalen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ca, bdxtail, bxtca);
                int temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(bxtcalen, bxtca, 2.0 * bdx, temp16a);
                Span<double> bxtaa = stackalloc double[8];
                int bxtaalen = ArithmeticFunctions.ScaleExpansionZeroelim(4, aa, bdxtail, bxtaa);
                int temp16blen = ArithmeticFunctions.ScaleExpansionZeroelim(bxtaalen, bxtaa, cdy, temp16b);
                Span<double> bxtcc = stackalloc double[8];
                int bxtcclen = ArithmeticFunctions.ScaleExpansionZeroelim(4, cc, bdxtail, bxtcc);
                int temp16clen = ArithmeticFunctions.ScaleExpansionZeroelim(bxtcclen, bxtcc, -ady, temp16c);
                int temp32alen = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                        temp16blen, temp16b, temp32a);
                int temp48len = ArithmeticFunctions.FastExpansionSumZeroelim(temp16clen, temp16c,
                                                        temp32alen, temp32a, temp48);
                finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                        temp48, finother);
                finswap = finnow;
                finnow = finother; 
                finother = finswap;
            }

            int bytcalen = 0;
            Span<double> bytca = stackalloc double[8];
            if (bdytail != 0.0)
            {
                bytcalen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ca, bdytail, bytca);
                int temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(bytcalen, bytca, 2.0 * bdy,
                                                      temp16a);

                Span<double> bytcc = stackalloc double[8];
                int bytcclen = ArithmeticFunctions.ScaleExpansionZeroelim(4, cc, bdytail, bytcc);
                int temp16blen = ArithmeticFunctions.ScaleExpansionZeroelim(bytcclen, bytcc, adx, temp16b);

                Span<double> bytaa = stackalloc double[8];
                int bytaalen = ArithmeticFunctions.ScaleExpansionZeroelim(4, aa, bdytail, bytaa);
                int temp16clen = ArithmeticFunctions.ScaleExpansionZeroelim(bytaalen, bytaa, -cdx, temp16c);
                int temp32alen = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                        temp16blen, temp16b, temp32a);
                int temp48len = ArithmeticFunctions.FastExpansionSumZeroelim(temp16clen, temp16c,
                                                        temp32alen, temp32a, temp48);
                finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                        temp48, finother);
                finswap = finnow;
                finnow = finother; 
                finother = finswap;
            }

            int cxtablen = 0;
            Span<double> cxtab = stackalloc double[8];
            if (cdxtail != 0.0)
            {
                cxtablen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ab, cdxtail, cxtab);
                var temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(cxtablen, cxtab, 2.0 * cdx,
                                                      temp16a);
                Span<double> cxtbb = stackalloc double[8];
                int cxtbblen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bb, cdxtail, cxtbb);
                int temp16blen = ArithmeticFunctions.ScaleExpansionZeroelim(cxtbblen, cxtbb, ady, temp16b);
                Span<double> cxtaa = stackalloc double[8];
                int cxtaalen = ArithmeticFunctions.ScaleExpansionZeroelim(4, aa, cdxtail, cxtaa);
                int temp16clen = ArithmeticFunctions.ScaleExpansionZeroelim(cxtaalen, cxtaa, -bdy, temp16c);

                int temp32alen = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                        temp16blen, temp16b, temp32a);
                int temp48len = ArithmeticFunctions.FastExpansionSumZeroelim(temp16clen, temp16c,
                                                        temp32alen, temp32a, temp48);
                finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                        temp48, finother);
                finswap = finnow;
                finnow = finother; 
                finother = finswap;
            }

            int cytablen = 0;
            Span<double> cytab = stackalloc double[8];
            if (cdytail != 0.0)
            {
                cytablen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ab, cdytail, cytab);
                int temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(cytablen, cytab, 2.0 * cdy,
                                                      temp16a);

                Span<double> cytaa = stackalloc double[8];
                int cytaalen = ArithmeticFunctions.ScaleExpansionZeroelim(4, aa, cdytail, cytaa);
                int temp16blen = ArithmeticFunctions.ScaleExpansionZeroelim(cytaalen, cytaa, bdx, temp16b);

                Span<double> cytbb = stackalloc double[8];
                int cytbblen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bb, cdytail, cytbb);
                int temp16clen = ArithmeticFunctions.ScaleExpansionZeroelim(cytbblen, cytbb, -adx, temp16c);

                int temp32alen = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                        temp16blen, temp16b, temp32a);
                int temp48len = ArithmeticFunctions.FastExpansionSumZeroelim(temp16clen, temp16c,
                                                        temp32alen, temp32a, temp48);
                finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                        temp48, finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap;
            }

            int bcttlen;
            Span<double> bct = stackalloc double[8];
            int bctlen;
            Span<double> bctt = stackalloc double[4];
            Span<double> temp64 = stackalloc double[64];
            Span<double> u = stackalloc double[5];
            Span<double> v = stackalloc double[5];
            if ((adxtail != 0.0) || (adytail != 0.0))
            {
                if ((bdxtail != 0.0) || (bdytail != 0.0)
                    || (cdxtail != 0.0) || (cdytail != 0.0))
                {
                    MacrosHelpers.TwoProduct(bdxtail, cdy, out double ti1, out double ti0);
                    MacrosHelpers.TwoProduct(bdx, cdytail, out double tj1, out double tj0);
                    MacrosHelpers.TwoTwoSum(ti1, ti0, tj1, tj0, out u[3], out u[2], out u[1], out u[0]);
                    double negate = -bdy;
                    MacrosHelpers.TwoProduct(cdxtail, negate, out ti1, out ti0);
                    negate = -bdytail;
                    MacrosHelpers.TwoProduct(cdx, negate, out tj1, out tj0);

                    MacrosHelpers.TwoTwoSum(ti1, ti0, tj1, tj0, out v[3], out v[2], out v[1], out v[0]);
                    bctlen = ArithmeticFunctions.FastExpansionSumZeroelim(4, u, 4, v, bct);

                    MacrosHelpers.TwoProduct(bdxtail, cdytail, out ti1, out ti0);
                    MacrosHelpers.TwoProduct(cdxtail, bdytail, out tj1, out tj0);
                    MacrosHelpers.TwoTwoDiff(ti1, ti0, tj1, tj0, out bctt[3], out bctt[2], out bctt[1], out bctt[0]);
                    bcttlen = 4;
                }
                else
                {
                    bct[0] = 0.0;
                    bctlen = 1;
                    bctt[0] = 0.0;
                    bcttlen = 1;
                }

                Span<double> temp8 = stackalloc double[8];
                if (adxtail != 0.0)
                {
                    int temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(axtbclen, axtbc, adxtail, temp16a);
                    Span<double> axtbct = stackalloc double[16];
                    int axtbctlen = ArithmeticFunctions.ScaleExpansionZeroelim(bctlen, bct, adxtail, axtbct);
                    int temp32alen = ArithmeticFunctions.ScaleExpansionZeroelim(axtbctlen, axtbct, 2.0 * adx,
                                                          temp32a);
                    int temp48len = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp32alen, temp32a, temp48);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                            temp48, finother);
                    finswap = finnow; 
                    finnow = finother;
                    finother = finswap;
                    if (bdytail != 0.0)
                    {
                        int temp8len = ArithmeticFunctions.ScaleExpansionZeroelim(4, cc, adxtail, temp8);
                        temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(temp8len, temp8, bdytail,
                                                              temp16a);
                        finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp16alen,
                                                                temp16a, finother);
                        finswap = finnow;
                        finnow = finother; 
                        finother = finswap;
                    }
                    if (cdytail != 0.0)
                    {
                        int temp8len = ArithmeticFunctions.ScaleExpansionZeroelim(4, bb, -adxtail, temp8);
                        temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(temp8len, temp8, cdytail,
                                                              temp16a);
                        finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp16alen,
                                                                temp16a, finother);
                        finswap = finnow;
                        finnow = finother;
                        finother = finswap;
                    }

                    temp32alen = ArithmeticFunctions.ScaleExpansionZeroelim(axtbctlen, axtbct, adxtail,
                                                          temp32a);
                    Span<double> axtbctt = stackalloc double[8];
                    int axtbcttlen = ArithmeticFunctions.ScaleExpansionZeroelim(bcttlen, bctt, adxtail, axtbctt);
                    temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(axtbcttlen, axtbctt, 2.0 * adx,
                                                          temp16a);
                    int temp16blen = ArithmeticFunctions.ScaleExpansionZeroelim(axtbcttlen, axtbctt, adxtail,
                                                          temp16b);
                    int temp32blen = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp16blen, temp16b, temp32b);
                    int temp64len = ArithmeticFunctions.FastExpansionSumZeroelim(temp32alen, temp32a,
                                                            temp32blen, temp32b, temp64);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp64len,
                                                            temp64, finother);
                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;
                }
                if (adytail != 0.0)
                {
                    Span<double> aytbct = stackalloc double[16];
                    int temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(aytbclen, aytbc, adytail, temp16a);
                    int aytbctlen = ArithmeticFunctions.ScaleExpansionZeroelim(bctlen, bct, adytail, aytbct);
                    int temp32alen = ArithmeticFunctions.ScaleExpansionZeroelim(aytbctlen, aytbct, 2.0 * ady,
                                                          temp32a);
                    int temp48len = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp32alen, temp32a, temp48);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                            temp48, finother);

                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;

                    temp32alen = ArithmeticFunctions.ScaleExpansionZeroelim(aytbctlen, aytbct, adytail,
                                                          temp32a);
                    Span<double> aytbctt = stackalloc double[16];
                    int aytbcttlen = ArithmeticFunctions.ScaleExpansionZeroelim(bcttlen, bctt, adytail, aytbctt);
                    temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(aytbcttlen, aytbctt, 2.0 * ady,
                                                          temp16a);
                    int temp16blen = ArithmeticFunctions.ScaleExpansionZeroelim(aytbcttlen, aytbctt, adytail,
                                                          temp16b);
                    int temp32blen = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp16blen, temp16b, temp32b);
                    int temp64len = ArithmeticFunctions.FastExpansionSumZeroelim(temp32alen, temp32a,
                                                            temp32blen, temp32b, temp64);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp64len,
                                                            temp64, finother);
                    finswap = finnow; 
                    finnow = finother;
                    finother = finswap;
                }
            }
            Span<double> cat = stackalloc double[8];
            Span<double> catt = stackalloc double[4];
            if ((bdxtail != 0.0) || (bdytail != 0.0))
            {
                int cattlen;
                int catlen;
                if ((cdxtail != 0.0) || (cdytail != 0.0)
                    || (adxtail != 0.0) || (adytail != 0.0))
                {
                    MacrosHelpers.TwoProduct(cdxtail, ady, out double ti1, out double ti0);
                    MacrosHelpers.TwoProduct(cdx, adytail, out double tj1, out double tj0);

                    MacrosHelpers.TwoTwoSum(ti1, ti0, tj1, tj0, out u[3], out u[2], out u[1], out u[0]);

                    double negate = -cdy;
                    MacrosHelpers.TwoProduct(adxtail, negate, out ti1, out ti0);
                    negate = -cdytail;
                    MacrosHelpers.TwoProduct(adx, negate, out tj1, out tj0);
                    MacrosHelpers.TwoTwoSum(ti1, ti0, tj1, tj0, out v[3], out v[2], out v[1], out v[0]);

                    catlen = ArithmeticFunctions.FastExpansionSumZeroelim(4, u, 4, v, cat);

                    MacrosHelpers.TwoProduct(cdxtail, adytail, out ti1, out ti0);
                    MacrosHelpers.TwoProduct(adxtail, cdytail, out tj1, out tj0);
                    MacrosHelpers.TwoTwoDiff(ti1, ti0, tj1, tj0, out catt[3], out catt[2], out catt[1], out catt[0]);
                    cattlen = 4;
                }
                else
                {
                    cat[0] = 0.0;
                    catlen = 1;
                    catt[0] = 0.0;
                    cattlen = 1;
                }

                Span<double> temp8 = stackalloc double[8];
                Span<double> bxtcat = stackalloc double[16];
                Span<double> bxtcatt = stackalloc double[8];
                if (bdxtail != 0.0)
                {
                    int temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(bxtcalen, bxtca, bdxtail, temp16a);
                    int bxtcatlen = ArithmeticFunctions.ScaleExpansionZeroelim(catlen, cat, bdxtail, bxtcat);
                    int temp32alen = ArithmeticFunctions.ScaleExpansionZeroelim(bxtcatlen, bxtcat, 2.0 * bdx,
                                                          temp32a);
                    int temp48len = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                             temp32alen, temp32a, temp48);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                            temp48, finother);
                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;
                    if (cdytail != 0.0)
                    {
                        int temp8len = ArithmeticFunctions.ScaleExpansionZeroelim(4, aa, bdxtail, temp8);
                        temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(temp8len, temp8, cdytail,
                                                              temp16a);
                        finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp16alen,
                                                                temp16a, finother);
                        finswap = finnow;
                        finnow = finother;
                        finother = finswap;
                    }
                    if (adytail != 0.0)
                    {
                        int temp8len = ArithmeticFunctions.ScaleExpansionZeroelim(4, cc, -bdxtail, temp8);
                        temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(temp8len, temp8, adytail,
                                                              temp16a);
                        finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp16alen,
                                                                temp16a, finother);
                        finswap = finnow; 
                        finnow = finother;
                        finother = finswap;
                    }

                    temp32alen = ArithmeticFunctions.ScaleExpansionZeroelim(bxtcatlen, bxtcat, bdxtail,
                                                          temp32a);
                    int bxtcattlen = ArithmeticFunctions.ScaleExpansionZeroelim(cattlen, catt, bdxtail, bxtcatt);
                    temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(bxtcattlen, bxtcatt, 2.0 * bdx,
                                                          temp16a);
                    int temp16blen = ArithmeticFunctions.ScaleExpansionZeroelim(bxtcattlen, bxtcatt, bdxtail,
                                                          temp16b);
                    int temp32blen = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp16blen, temp16b, temp32b);
                    int temp64len = ArithmeticFunctions.FastExpansionSumZeroelim(temp32alen, temp32a,
                                                            temp32blen, temp32b, temp64);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp64len,
                                                            temp64, finother);
                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;
                }
                Span<double> bytcat = stackalloc double[16];
                Span<double> bytcatt = stackalloc double[8];
                if (bdytail != 0.0)
                {
                    int temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(bytcalen, bytca, bdytail, temp16a);
                    int bytcatlen = ArithmeticFunctions.ScaleExpansionZeroelim(catlen, cat, bdytail, bytcat);
                    int temp32alen = ArithmeticFunctions.ScaleExpansionZeroelim(bytcatlen, bytcat, 2.0 * bdy,
                                                          temp32a);
                    int temp48len = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp32alen, temp32a, temp48);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                            temp48, finother);
                    finswap = finnow;
                    finnow = finother; 
                    finother = finswap;


                    temp32alen = ArithmeticFunctions.ScaleExpansionZeroelim(bytcatlen, bytcat, bdytail,
                                                          temp32a);
                    int bytcattlen = ArithmeticFunctions.ScaleExpansionZeroelim(cattlen, catt, bdytail, bytcatt);
                    temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(bytcattlen, bytcatt, 2.0 * bdy,
                                                          temp16a);
                    int temp16blen = ArithmeticFunctions.ScaleExpansionZeroelim(bytcattlen, bytcatt, bdytail,
                                                          temp16b);
                    int temp32blen = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp16blen, temp16b, temp32b);
                    int temp64len = ArithmeticFunctions.FastExpansionSumZeroelim(temp32alen, temp32a,
                                                            temp32blen, temp32b, temp64);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp64len,
                                                            temp64, finother);
                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;
                }
            }
            Span<double> abt = stackalloc double[8];
            Span<double> abtt = stackalloc double[4];
            if ((cdxtail != 0.0) || (cdytail != 0.0))
            {
                int abtlen;
                int abttlen;
                if ((adxtail != 0.0) || (adytail != 0.0)
                    || (bdxtail != 0.0) || (bdytail != 0.0))
                {
                    MacrosHelpers.TwoProduct(adxtail, bdy, out double ti1, out double ti0);
                    MacrosHelpers.TwoProduct(adx, bdytail, out double tj1, out double tj0);
                    MacrosHelpers.TwoTwoSum(ti1, ti0, tj1, tj0, out u[3], out u[2], out u[1], out u[0]);
                    double negate = -ady;
                    MacrosHelpers.TwoProduct(bdxtail, negate, out ti1, out ti0);
                    negate = -adytail;
                    MacrosHelpers.TwoProduct(bdx, negate, out tj1, out tj0);
                    MacrosHelpers.TwoTwoSum(ti1, ti0, tj1, tj0, out v[3], out v[2], out v[1], out v[0]);
                    abtlen = ArithmeticFunctions.FastExpansionSumZeroelim(4, u, 4, v, abt);

                    MacrosHelpers.TwoProduct(adxtail, bdytail, out ti1, out ti0);
                    MacrosHelpers.TwoProduct(bdxtail, adytail, out tj1, out tj0);
                    MacrosHelpers.TwoTwoDiff(ti1, ti0, tj1, tj0, out abtt[3], out abtt[2], out abtt[1], out abtt[0]);
                    abttlen = 4;
                }
                else
                {
                    abt[0] = 0.0;
                    abtlen = 1;
                    abtt[0] = 0.0;
                    abttlen = 1;
                }

                Span<double> temp8 = stackalloc double[8];
                Span<double> cxtabt = stackalloc double[16];
                Span<double> cxtabtt = stackalloc double[8];
                if (cdxtail != 0.0)
                {
                    int temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(cxtablen, cxtab, cdxtail, temp16a);
                    int cxtabtlen = ArithmeticFunctions.ScaleExpansionZeroelim(abtlen, abt, cdxtail, cxtabt);
                    int temp32alen = ArithmeticFunctions.ScaleExpansionZeroelim(cxtabtlen, cxtabt, 2.0 * cdx,
                                                          temp32a);
                    int temp48len = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp32alen, temp32a, temp48);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                            temp48, finother);
                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;

                    if (adytail != 0.0)
                    {
                        int temp8len = ArithmeticFunctions.ScaleExpansionZeroelim(4, bb, cdxtail, temp8);
                        temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(temp8len, temp8, adytail,
                                                              temp16a);
                        finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp16alen,
                                                                temp16a, finother);
                        finswap = finnow; 
                        finnow = finother; 
                        finother = finswap;
                    }
                    if (bdytail != 0.0)
                    {
                        int temp8len = ArithmeticFunctions.ScaleExpansionZeroelim(4, aa, -cdxtail, temp8);
                        temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(temp8len, temp8, bdytail,
                                                              temp16a);
                        finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp16alen,
                                                                temp16a, finother);
                        finswap = finnow; 
                        finnow = finother; 
                        finother = finswap;
                    }

                    temp32alen = ArithmeticFunctions.ScaleExpansionZeroelim(cxtabtlen, cxtabt, cdxtail,
                                                          temp32a);
                    int cxtabttlen = ArithmeticFunctions.ScaleExpansionZeroelim(abttlen, abtt, cdxtail, cxtabtt);
                    temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(cxtabttlen, cxtabtt, 2.0 * cdx,
                                                          temp16a);
                    int temp16blen = ArithmeticFunctions.ScaleExpansionZeroelim(cxtabttlen, cxtabtt, cdxtail,
                                                          temp16b);
                    int temp32blen = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp16blen, temp16b, temp32b);
                    int temp64len = ArithmeticFunctions.FastExpansionSumZeroelim(temp32alen, temp32a,
                                                            temp32blen, temp32b, temp64);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp64len,
                                                            temp64, finother);
                    finswap = finnow; 
                    finnow = finother; 
                    finother = finswap;
                }

                Span<double> cytabt = stackalloc double[16];
                Span<double> cytabtt = stackalloc double[8];
                if (cdytail != 0.0)
                {
                    int temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(cytablen, cytab, cdytail, temp16a);
                    int cytabtlen = ArithmeticFunctions.ScaleExpansionZeroelim(abtlen, abt, cdytail, cytabt);
                    int temp32alen = ArithmeticFunctions.ScaleExpansionZeroelim(cytabtlen, cytabt, 2.0 * cdy,
                                                          temp32a);
                    int temp48len = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp32alen, temp32a, temp48);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                            temp48, finother);
                    finswap = finnow; finnow = finother; finother = finswap;


                    temp32alen = ArithmeticFunctions.ScaleExpansionZeroelim(cytabtlen, cytabt, cdytail,
                                                          temp32a);
                    int cytabttlen = ArithmeticFunctions.ScaleExpansionZeroelim(abttlen, abtt, cdytail, cytabtt);
                    temp16alen = ArithmeticFunctions.ScaleExpansionZeroelim(cytabttlen, cytabtt, 2.0 * cdy,
                                                          temp16a);
                    int temp16blen = ArithmeticFunctions.ScaleExpansionZeroelim(cytabttlen, cytabtt, cdytail,
                                                          temp16b);
                    int temp32blen = ArithmeticFunctions.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp16blen, temp16b, temp32b);
                    int temp64len = ArithmeticFunctions.FastExpansionSumZeroelim(temp32alen, temp32a,
                                                            temp32blen, temp32b, temp64);
                    finlength = ArithmeticFunctions.FastExpansionSumZeroelim(finlength, finnow, temp64len,
                                                            temp64, finother);
                    finswap = finnow; 
                    finnow = finother; 
                    finother = finswap;
                }
            }

            return finnow[finlength - 1];
        }

        public static double Robust(ReadOnlySpan<double> pa, ReadOnlySpan<double> pb, ReadOnlySpan<double> pc, ReadOnlySpan<double> pd)
        {
            double adx = pa[0] - pd[0];
            double bdx = pb[0] - pd[0];
            double cdx = pc[0] - pd[0];
            double ady = pa[1] - pd[1];
            double bdy = pb[1] - pd[1];
            double cdy = pc[1] - pd[1];

            double bdxcdy = bdx * cdy;
            double cdxbdy = cdx * bdy;
            double alift = adx * adx + ady * ady;

            double cdxady = cdx * ady;
            double adxcdy = adx * cdy;
            double blift = bdx * bdx + bdy * bdy;

            double adxbdy = adx * bdy;
            double bdxady = bdx * ady;
            double clift = cdx * cdx + cdy * cdy;

            double det = alift * (bdxcdy - cdxbdy)
                  + blift * (cdxady - adxcdy)
                  + clift * (adxbdy - bdxady);

            double permanent = (Math.Abs(bdxcdy) + Math.Abs(cdxbdy)) * alift
                        + (Math.Abs(cdxady) + Math.Abs(adxcdy)) * blift
                        + (Math.Abs(adxbdy) + Math.Abs(bdxady)) * clift;

            double errbound = MacrosHelpers.IccerrboundA * permanent;
            if ((det > errbound) || (-det > errbound))
            {
                return det;
            }

            return Adapt(pa, pb, pc, pd, permanent);
        }
    }
}
