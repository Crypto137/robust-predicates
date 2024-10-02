using System;

namespace RobustPredicates
{
    public static class InSphere
    {
        public static unsafe double Fast(double* pa, double* pb, double* pc, double* pd, double* pe)
        {
            double aex = pa[0] - pe[0];
            double bex = pb[0] - pe[0];
            double cex = pc[0] - pe[0];
            double dex = pd[0] - pe[0];
            double aey = pa[1] - pe[1];
            double bey = pb[1] - pe[1];
            double cey = pc[1] - pe[1];
            double dey = pd[1] - pe[1];
            double aez = pa[2] - pe[2];
            double bez = pb[2] - pe[2];
            double cez = pc[2] - pe[2];
            double dez = pd[2] - pe[2];

            double ab = aex * bey - bex * aey;
            double bc = bex * cey - cex * bey;
            double cd = cex * dey - dex * cey;
            double da = dex * aey - aex * dey;

            double ac = aex * cey - cex * aey;
            double bd = bex * dey - dex * bey;

            double abc = aez * bc - bez * ac + cez * ab;
            double bcd = bez * cd - cez * bd + dez * bc;
            double cda = cez * da + dez * ac + aez * cd;
            double dab = dez * ab + aez * bd + bez * da;

            double alift = aex * aex + aey * aey + aez * aez;
            double blift = bex * bex + bey * bey + bez * bez;
            double clift = cex * cex + cey * cey + cez * cez;
            double dlift = dex * dex + dey * dey + dez * dez;

            return (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
        }

        public static unsafe double Fast(ReadOnlySpan<double> pa, ReadOnlySpan<double> pb, ReadOnlySpan<double> pc, ReadOnlySpan<double> pd, ReadOnlySpan<double> pe)
        {
            fixed (double* paPtr = pa)
            fixed (double* pbPtr = pb)
            fixed (double* pcPtr = pc)
            fixed (double* pdPtr = pd)
            fixed (double* pePtr = pe)
            {
                return Fast(paPtr, pbPtr, pcPtr, pdPtr, pePtr);
            }
        }

        public static unsafe double Robust(double* pa, double* pb, double* pc, double* pd, double* pe)
        {
            double aex, bex, cex, dex;
            double aey, bey, cey, dey;
            double aez, bez, cez, dez;
            double aexbey, bexaey, bexcey, cexbey, cexdey, dexcey, dexaey, aexdey;
            double aexcey, cexaey, bexdey, dexbey;
            double alift, blift, clift, dlift;
            double ab, bc, cd, da, ac, bd;
            double abc, bcd, cda, dab;
            double aezplus, bezplus, cezplus, dezplus;
            double aexbeyplus, bexaeyplus, bexceyplus, cexbeyplus;
            double cexdeyplus, dexceyplus, dexaeyplus, aexdeyplus;
            double aexceyplus, cexaeyplus, bexdeyplus, dexbeyplus;
            double det;
            double permanent, errbound;

            aex = pa[0] - pe[0];
            bex = pb[0] - pe[0];
            cex = pc[0] - pe[0];
            dex = pd[0] - pe[0];
            aey = pa[1] - pe[1];
            bey = pb[1] - pe[1];
            cey = pc[1] - pe[1];
            dey = pd[1] - pe[1];
            aez = pa[2] - pe[2];
            bez = pb[2] - pe[2];
            cez = pc[2] - pe[2];
            dez = pd[2] - pe[2];

            aexbey = aex * bey;
            bexaey = bex * aey;
            ab = aexbey - bexaey;
            bexcey = bex * cey;
            cexbey = cex * bey;
            bc = bexcey - cexbey;
            cexdey = cex * dey;
            dexcey = dex * cey;
            cd = cexdey - dexcey;
            dexaey = dex * aey;
            aexdey = aex * dey;
            da = dexaey - aexdey;

            aexcey = aex * cey;
            cexaey = cex * aey;
            ac = aexcey - cexaey;
            bexdey = bex * dey;
            dexbey = dex * bey;
            bd = bexdey - dexbey;

            abc = aez * bc - bez * ac + cez * ab;
            bcd = bez * cd - cez * bd + dez * bc;
            cda = cez * da + dez * ac + aez * cd;
            dab = dez * ab + aez * bd + bez * da;

            alift = aex * aex + aey * aey + aez * aez;
            blift = bex * bex + bey * bey + bez * bez;
            clift = cex * cex + cey * cey + cez * cez;
            dlift = dex * dex + dey * dey + dez * dez;

            det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

            aezplus = Math.Abs(aez);
            bezplus = Math.Abs(bez);
            cezplus = Math.Abs(cez);
            dezplus = Math.Abs(dez);
            aexbeyplus = Math.Abs(aexbey);
            bexaeyplus = Math.Abs(bexaey);
            bexceyplus = Math.Abs(bexcey);
            cexbeyplus = Math.Abs(cexbey);
            cexdeyplus = Math.Abs(cexdey);
            dexceyplus = Math.Abs(dexcey);
            dexaeyplus = Math.Abs(dexaey);
            aexdeyplus = Math.Abs(aexdey);
            aexceyplus = Math.Abs(aexcey);
            cexaeyplus = Math.Abs(cexaey);
            bexdeyplus = Math.Abs(bexdey);
            dexbeyplus = Math.Abs(dexbey);
            permanent = ((cexdeyplus + dexceyplus) * bezplus
                         + (dexbeyplus + bexdeyplus) * cezplus
                         + (bexceyplus + cexbeyplus) * dezplus)
                      * alift
                      + ((dexaeyplus + aexdeyplus) * cezplus
                         + (aexceyplus + cexaeyplus) * dezplus
                         + (cexdeyplus + dexceyplus) * aezplus)
                      * blift
                      + ((aexbeyplus + bexaeyplus) * dezplus
                         + (bexdeyplus + dexbeyplus) * aezplus
                         + (dexaeyplus + aexdeyplus) * bezplus)
                      * clift
                      + ((bexceyplus + cexbeyplus) * aezplus
                         + (cexaeyplus + aexceyplus) * bezplus
                         + (aexbeyplus + bexaeyplus) * cezplus)
                      * dlift;

            errbound = MacrosHelpers.IsperrboundA * permanent;
            if ((det > errbound) || (-det > errbound))
            {
                return det;
            }

            return Adapt(pa, pb, pc, pd, pe, permanent);
        }

        public static unsafe double Robust(ReadOnlySpan<double> pa, ReadOnlySpan<double> pb, ReadOnlySpan<double> pc, ReadOnlySpan<double> pd, ReadOnlySpan<double> pe)
        {
            fixed (double* paPtr = pa)
            fixed (double* pbPtr = pb)
            fixed (double* pcPtr = pc)
            fixed (double* pdPtr = pd)
            fixed (double* pePtr = pe)
            {
                return Robust(paPtr, pbPtr, pcPtr, pdPtr, pePtr);
            }
        }

        private static unsafe double Adapt(double* pa, double* pb, double* pc, double* pd, double* pe, double permanent)
        {
            double aex = pa[0] - pe[0];
            double bex = pb[0] - pe[0];
            double cex = pc[0] - pe[0];
            double dex = pd[0] - pe[0];
            double aey = pa[1] - pe[1];
            double bey = pb[1] - pe[1];
            double cey = pc[1] - pe[1];
            double dey = pd[1] - pe[1];
            double aez = pa[2] - pe[2];
            double bez = pb[2] - pe[2];
            double cez = pc[2] - pe[2];
            double dez = pd[2] - pe[2];

            MacrosHelpers.TwoProduct(aex, bey, out double aexbey1, out double aexbey0);
            MacrosHelpers.TwoProduct(bex, aey, out double bexaey1, out double bexaey0);
            double* ab = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(aexbey1, aexbey0, bexaey1, bexaey0, out ab[3], out ab[2], out ab[1], out ab[0]);

            MacrosHelpers.TwoProduct(bex, cey, out double bexcey1, out double bexcey0);
            MacrosHelpers.TwoProduct(cex, bey, out double cexbey1, out double cexbey0);
            double* bc = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(bexcey1, bexcey0, cexbey1, cexbey0, out bc[3], out bc[2], out bc[1], out bc[0]);

            MacrosHelpers.TwoProduct(cex, dey, out double cexdey1, out double cexdey0);
            MacrosHelpers.TwoProduct(dex, cey, out double dexcey1, out double dexcey0);
            double* cd = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(cexdey1, cexdey0, dexcey1, dexcey0, out cd[3], out cd[2], out cd[1], out cd[0]);

            MacrosHelpers.TwoProduct(dex, aey, out double dexaey1, out double dexaey0);
            MacrosHelpers.TwoProduct(aex, dey, out double aexdey1, out double aexdey0);
            double* da = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(dexaey1, dexaey0, aexdey1, aexdey0, out da[3], out da[2], out da[1], out da[0]);

            MacrosHelpers.TwoProduct(aex, cey, out double aexcey1, out double aexcey0);
            MacrosHelpers.TwoProduct(cex, aey, out double cexaey1, out double cexaey0);
            double* ac = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(aexcey1, aexcey0, cexaey1, cexaey0, out ac[3], out ac[2], out ac[1], out ac[0]);

            MacrosHelpers.TwoProduct(bex, dey, out double bexdey1, out double bexdey0);
            MacrosHelpers.TwoProduct(dex, bey, out double dexbey1, out double dexbey0);
            double* bd = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(bexdey1, bexdey0, dexbey1, dexbey0, out bd[3], out bd[2], out bd[1], out bd[0]);

            double* temp8a = stackalloc double[8];
            double* temp8b = stackalloc double[8];
            double* temp8c = stackalloc double[8];
            int temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, cd, bez, temp8a);
            int temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bd, -cez, temp8b);
            int temp8clen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bc, dez, temp8c);
            double* temp16 = stackalloc double[16];
            int temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a,
                                                    temp8blen, temp8b, temp16);
            double* temp24 = stackalloc double[24];
            int temp24len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8clen, temp8c,
                                                    temp16len, temp16, temp24);

            double* temp48 = stackalloc double[48];
            int temp48len = ArithmeticFunctions.ScaleExpansionZeroelim(temp24len, temp24, aex, temp48);
            double* xdet = stackalloc double[96];
            int xlen = ArithmeticFunctions.ScaleExpansionZeroelim(temp48len, temp48, -aex, xdet);
            temp48len = ArithmeticFunctions.ScaleExpansionZeroelim(temp24len, temp24, aey, temp48);
            double* ydet = stackalloc double[96];
            int ylen = ArithmeticFunctions.ScaleExpansionZeroelim(temp48len, temp48, -aey, ydet);
            temp48len = ArithmeticFunctions.ScaleExpansionZeroelim(temp24len, temp24, aez, temp48);
            double* zdet = stackalloc double[96];
            var zlen = ArithmeticFunctions.ScaleExpansionZeroelim(temp48len, temp48, -aez, zdet);

            double* xydet = stackalloc double[192];
            int xylen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, xdet, ylen, ydet, xydet);
            double* adet = stackalloc double[288];
            int alen = ArithmeticFunctions.FastExpansionSumZeroelim(xylen, xydet, zlen, zdet, adet);

            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, da, cez, temp8a);
            temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ac, dez, temp8b);
            temp8clen = ArithmeticFunctions.ScaleExpansionZeroelim(4, cd, aez, temp8c);
            temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a,
                                                    temp8blen, temp8b, temp16);
            temp24len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8clen, temp8c,
                                                    temp16len, temp16, temp24);
            temp48len = ArithmeticFunctions.ScaleExpansionZeroelim(temp24len, temp24, bex, temp48);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(temp48len, temp48, bex, xdet);
            temp48len = ArithmeticFunctions.ScaleExpansionZeroelim(temp24len, temp24, bey, temp48);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(temp48len, temp48, bey, ydet);
            temp48len = ArithmeticFunctions.ScaleExpansionZeroelim(temp24len, temp24, bez, temp48);
            zlen = ArithmeticFunctions.ScaleExpansionZeroelim(temp48len, temp48, bez, zdet);
            xylen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, xdet, ylen, ydet, xydet);
            double* bdet = stackalloc double[288];
            int blen = ArithmeticFunctions.FastExpansionSumZeroelim(xylen, xydet, zlen, zdet, bdet);

            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ab, dez, temp8a);
            temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bd, aez, temp8b);
            temp8clen = ArithmeticFunctions.ScaleExpansionZeroelim(4, da, bez, temp8c);
            temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a,
                                                    temp8blen, temp8b, temp16);
            temp24len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8clen, temp8c,
                                                    temp16len, temp16, temp24);
            temp48len = ArithmeticFunctions.ScaleExpansionZeroelim(temp24len, temp24, cex, temp48);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(temp48len, temp48, -cex, xdet);
            temp48len = ArithmeticFunctions.ScaleExpansionZeroelim(temp24len, temp24, cey, temp48);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(temp48len, temp48, -cey, ydet);
            temp48len = ArithmeticFunctions.ScaleExpansionZeroelim(temp24len, temp24, cez, temp48);
            zlen = ArithmeticFunctions.ScaleExpansionZeroelim(temp48len, temp48, -cez, zdet);
            xylen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, xdet, ylen, ydet, xydet);
            double* cdet = stackalloc double[288];
            int clen = ArithmeticFunctions.FastExpansionSumZeroelim(xylen, xydet, zlen, zdet, cdet);

            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bc, aez, temp8a);
            temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ac, -bez, temp8b);
            temp8clen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ab, cez, temp8c);
            temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a,
                                                    temp8blen, temp8b, temp16);
            temp24len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8clen, temp8c,
                                                    temp16len, temp16, temp24);
            temp48len = ArithmeticFunctions.ScaleExpansionZeroelim(temp24len, temp24, dex, temp48);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(temp48len, temp48, dex, xdet);
            temp48len = ArithmeticFunctions.ScaleExpansionZeroelim(temp24len, temp24, dey, temp48);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(temp48len, temp48, dey, ydet);
            temp48len = ArithmeticFunctions.ScaleExpansionZeroelim(temp24len, temp24, dez, temp48);
            zlen = ArithmeticFunctions.ScaleExpansionZeroelim(temp48len, temp48, dez, zdet);
            xylen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, xdet, ylen, ydet, xydet);
            double* ddet = stackalloc double[288];
            int dlen = ArithmeticFunctions.FastExpansionSumZeroelim(xylen, xydet, zlen, zdet, ddet);

            double* abdet = stackalloc double[576];
            double* cddet = stackalloc double[576];
            double* fin1 = stackalloc double[1152];
            int ablen = ArithmeticFunctions.FastExpansionSumZeroelim(alen, adet, blen, bdet, abdet);
            int cdlen = ArithmeticFunctions.FastExpansionSumZeroelim(clen, cdet, dlen, ddet, cddet);
            int finlength = ArithmeticFunctions.FastExpansionSumZeroelim(ablen, abdet, cdlen, cddet, fin1);
            double det = ArithmeticFunctions.Estimate(finlength, fin1);
            double errbound = MacrosHelpers.IsperrboundB * permanent;
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            MacrosHelpers.TwoDiffTail(pa[0], pe[0], aex, out double aextail);
            MacrosHelpers.TwoDiffTail(pa[1], pe[1], aey, out double aeytail);
            MacrosHelpers.TwoDiffTail(pa[2], pe[2], aez, out double aeztail);
            MacrosHelpers.TwoDiffTail(pb[0], pe[0], bex, out double bextail);
            MacrosHelpers.TwoDiffTail(pb[1], pe[1], bey, out double beytail);
            MacrosHelpers.TwoDiffTail(pb[2], pe[2], bez, out double beztail);
            MacrosHelpers.TwoDiffTail(pc[0], pe[0], cex, out double cextail);
            MacrosHelpers.TwoDiffTail(pc[1], pe[1], cey, out double ceytail);
            MacrosHelpers.TwoDiffTail(pc[2], pe[2], cez, out double ceztail);
            MacrosHelpers.TwoDiffTail(pd[0], pe[0], dex, out double dextail);
            MacrosHelpers.TwoDiffTail(pd[1], pe[1], dey, out double deytail);
            MacrosHelpers.TwoDiffTail(pd[2], pe[2], dez, out double deztail);

            if ((aextail == 0.0) && (aeytail == 0.0) && (aeztail == 0.0)
                && (bextail == 0.0) && (beytail == 0.0) && (beztail == 0.0)
                && (cextail == 0.0) && (ceytail == 0.0) && (ceztail == 0.0)
                && (dextail == 0.0) && (deytail == 0.0) && (deztail == 0.0))
            {
                return det;
            }

            double da3 = da[3];
            double ac3 = ac[3];
            double cd3 = cd[3];
            double ab3 = ab[3];
            double bc3 = bc[3];
            double bd3 = bd[3];

            errbound = MacrosHelpers.IsperrboundC * permanent + MacrosHelpers.Resulterrbound * Math.Abs(det);
            double abeps = (aex * beytail + bey * aextail)
                  - (aey * bextail + bex * aeytail);
            double bceps = (bex * ceytail + cey * bextail)
                  - (bey * cextail + cex * beytail);
            double cdeps = (cex * deytail + dey * cextail)
                  - (cey * dextail + dex * ceytail);
            double daeps = (dex * aeytail + aey * dextail)
                  - (dey * aextail + aex * deytail);
            double aceps = (aex * ceytail + cey * aextail)
                  - (aey * cextail + cex * aeytail);
            double bdeps = (bex * deytail + dey * bextail)
                  - (bey * dextail + dex * beytail);
            det += (((bex * bex + bey * bey + bez * bez)
                     * ((cez * daeps + dez * aceps + aez * cdeps)
                        + (ceztail * da3 + deztail * ac3 + aeztail * cd3))
                     + (dex * dex + dey * dey + dez * dez)
                     * ((aez * bceps - bez * aceps + cez * abeps)
                        + (aeztail * bc3 - beztail * ac3 + ceztail * ab3)))
                    - ((aex * aex + aey * aey + aez * aez)
                     * ((bez * cdeps - cez * bdeps + dez * bceps)
                        + (beztail * cd3 - ceztail * bd3 + deztail * bc3))
                     + (cex * cex + cey * cey + cez * cez)
                     * ((dez * abeps + aez * bdeps + bez * daeps)
                        + (deztail * ab3 + aeztail * bd3 + beztail * da3))))
                 + 2.0 * (((bex * bextail + bey * beytail + bez * beztail)
                           * (cez * da3 + dez * ac3 + aez * cd3)
                           + (dex * dextail + dey * deytail + dez * deztail)
                           * (aez * bc3 - bez * ac3 + cez * ab3))
                          - ((aex * aextail + aey * aeytail + aez * aeztail)
                           * (bez * cd3 - cez * bd3 + dez * bc3)
                           + (cex * cextail + cey * ceytail + cez * ceztail)
                           * (dez * ab3 + aez * bd3 + bez * da3)));
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            return Exact(pa, pb, pc, pd, pe);
        }

        private static unsafe double Exact(double* pa, double* pb, double* pc, double* pd, double* pe)
        {
            MacrosHelpers.TwoProduct(pa[0], pb[1], out double axby1, out double axby0);
            MacrosHelpers.TwoProduct(pb[0], pa[1], out double bxay1, out double bxay0);
            double* ab = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(axby1, axby0, bxay1, bxay0, out ab[3], out ab[2], out ab[1], out ab[0]);

            MacrosHelpers.TwoProduct(pb[0], pc[1], out double bxcy1, out double bxcy0);
            MacrosHelpers.TwoProduct(pc[0], pb[1], out double cxby1, out double cxby0);
            double* bc = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(bxcy1, bxcy0, cxby1, cxby0, out bc[3], out bc[2], out bc[1], out bc[0]);

            MacrosHelpers.TwoProduct(pc[0], pd[1], out double cxdy1, out double cxdy0);
            MacrosHelpers.TwoProduct(pd[0], pc[1], out double dxcy1, out double dxcy0);
            double* cd = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(cxdy1, cxdy0, dxcy1, dxcy0, out cd[3], out cd[2], out cd[1], out cd[0]);

            MacrosHelpers.TwoProduct(pd[0], pe[1], out double dxey1, out double dxey0);
            MacrosHelpers.TwoProduct(pe[0], pd[1], out double exdy1, out double exdy0);
            double* de = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(dxey1, dxey0, exdy1, exdy0, out de[3], out de[2], out de[1], out de[0]);

            MacrosHelpers.TwoProduct(pe[0], pa[1], out double exay1, out double exay0);
            MacrosHelpers.TwoProduct(pa[0], pe[1], out double axey1, out double axey0);
            double* ea = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(exay1, exay0, axey1, axey0, out ea[3], out ea[2], out ea[1], out ea[0]);

            MacrosHelpers.TwoProduct(pa[0], pc[1], out double axcy1, out double axcy0);
            MacrosHelpers.TwoProduct(pc[0], pa[1], out double cxay1, out double cxay0);
            double* ac = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(axcy1, axcy0, cxay1, cxay0, out ac[3], out ac[2], out ac[1], out ac[0]);

            MacrosHelpers.TwoProduct(pb[0], pd[1], out double bxdy1, out double bxdy0);
            MacrosHelpers.TwoProduct(pd[0], pb[1], out double dxby1, out double dxby0);
            double* bd = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(bxdy1, bxdy0, dxby1, dxby0, out bd[3], out bd[2], out bd[1], out bd[0]);

            MacrosHelpers.TwoProduct(pc[0], pe[1], out double cxey1, out double cxey0);
            MacrosHelpers.TwoProduct(pe[0], pc[1], out double excy1, out double excy0);
            double* ce = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(cxey1, cxey0, excy1, excy0, out ce[3], out ce[2], out ce[1], out ce[0]);

            MacrosHelpers.TwoProduct(pd[0], pa[1], out double dxay1, out double dxay0);
            MacrosHelpers.TwoProduct(pa[0], pd[1], out double axdy1, out double axdy0);
            double* da = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(dxay1, dxay0, axdy1, axdy0, out da[3], out da[2], out da[1], out da[0]);

            MacrosHelpers.TwoProduct(pe[0], pb[1], out double exby1, out double exby0);
            MacrosHelpers.TwoProduct(pb[0], pe[1], out double bxey1, out double bxey0);
            double* eb = stackalloc double[4];
            MacrosHelpers.TwoTwoDiff(exby1, exby0, bxey1, bxey0, out eb[3], out eb[2], out eb[1], out eb[0]);

            double* temp8a = stackalloc double[8];
            double* temp8b = stackalloc double[8];
            double* temp16 = stackalloc double[16];
            int temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bc, pa[2], temp8a);
            int temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ac, -pb[2], temp8b);
            int temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp8blen, temp8b,
                                                    temp16);
            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ab, pc[2], temp8a);
            double* abc = stackalloc double[24];
            int abclen = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp16len, temp16,
                                                 abc);

            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, cd, pb[2], temp8a);
            temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bd, -pc[2], temp8b);
            temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp8blen, temp8b,
                                                    temp16);
            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bc, pd[2], temp8a);
            double* bcd = stackalloc double[24];
            int bcdlen = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp16len, temp16,
                                                 bcd);

            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, de, pc[2], temp8a);
            temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ce, -pd[2], temp8b);
            temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp8blen, temp8b,
                                                    temp16);
            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, cd, pe[2], temp8a);
            double* cde = stackalloc double[24];
            int cdelen = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp16len, temp16,
                                                 cde);

            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ea, pd[2], temp8a);
            temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, da, -pe[2], temp8b);
            temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp8blen, temp8b,
                                                    temp16);
            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, de, pa[2], temp8a);
            double* dea = stackalloc double[24];
            int dealen = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp16len, temp16,
                                                 dea);

            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ab, pe[2], temp8a);
            temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, eb, -pa[2], temp8b);
            temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp8blen, temp8b,
                                                    temp16);
            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ea, pb[2], temp8a);
            double* eab = stackalloc double[24];
            int eablen = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp16len, temp16,
                                                 eab);

            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bd, pa[2], temp8a);
            temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, da, pb[2], temp8b);
            temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp8blen, temp8b,
                                                    temp16);
            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ab, pd[2], temp8a);
            double* abd = stackalloc double[24];
            int abdlen = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp16len, temp16,
                                                 abd);

            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ce, pb[2], temp8a);
            temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, eb, pc[2], temp8b);
            temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp8blen, temp8b,
                                                    temp16);
            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bc, pe[2], temp8a);
            double* bce = stackalloc double[24];
            int bcelen = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp16len, temp16,
                                                 bce);

            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, da, pc[2], temp8a);
            temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ac, pd[2], temp8b);
            temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp8blen, temp8b,
                                                    temp16);
            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, cd, pa[2], temp8a);
            double* cda = stackalloc double[24];
            int cdalen = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp16len, temp16,
                                                 cda);

            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, eb, pd[2], temp8a);
            temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, bd, pe[2], temp8b);
            temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp8blen, temp8b,
                                                    temp16);
            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, de, pb[2], temp8a);
            double* deb = stackalloc double[24];
            int deblen = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp16len, temp16,
                                                 deb);

            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ac, pe[2], temp8a);
            temp8blen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ce, pa[2], temp8b);
            temp16len = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp8blen, temp8b,
                                                    temp16);
            temp8alen = ArithmeticFunctions.ScaleExpansionZeroelim(4, ea, pc[2], temp8a);
            double* eac = stackalloc double[24];
            int eaclen = ArithmeticFunctions.FastExpansionSumZeroelim(temp8alen, temp8a, temp16len, temp16,
                                                 eac);

            double* temp48a = stackalloc double[48];
            double* temp48b = stackalloc double[48];
            int temp48alen = ArithmeticFunctions.FastExpansionSumZeroelim(cdelen, cde, bcelen, bce, temp48a);
            int temp48blen = ArithmeticFunctions.FastExpansionSumZeroelim(deblen, deb, bcdlen, bcd, temp48b);
            for (int i = 0; i < temp48blen; i++)
            {
                temp48b[i] = -temp48b[i];
            }

            double* bcde = stackalloc double[96];
            int bcdelen = ArithmeticFunctions.FastExpansionSumZeroelim(temp48alen, temp48a,
                                                  temp48blen, temp48b, bcde);

            double* temp192 = stackalloc double[192];
            int xlen = ArithmeticFunctions.ScaleExpansionZeroelim(bcdelen, bcde, pa[0], temp192);

            double* det384x = stackalloc double[384];
            double* det384y = stackalloc double[384];
            double* det384z = stackalloc double[384];
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(xlen, temp192, pa[0], det384x);
            int ylen = ArithmeticFunctions.ScaleExpansionZeroelim(bcdelen, bcde, pa[1], temp192);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(ylen, temp192, pa[1], det384y);
            int zlen = ArithmeticFunctions.ScaleExpansionZeroelim(bcdelen, bcde, pa[2], temp192);
            zlen = ArithmeticFunctions.ScaleExpansionZeroelim(zlen, temp192, pa[2], det384z);

            double* detxy = stackalloc double[786];
            double* adet = stackalloc double[1152];
            int xylen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, det384x, ylen, det384y, detxy);
            int alen = ArithmeticFunctions.FastExpansionSumZeroelim(xylen, detxy, zlen, det384z, adet);

            temp48alen = ArithmeticFunctions.FastExpansionSumZeroelim(dealen, dea, cdalen, cda, temp48a);
            temp48blen = ArithmeticFunctions.FastExpansionSumZeroelim(eaclen, eac, cdelen, cde, temp48b);
            for (int i = 0; i < temp48blen; i++)
            {
                temp48b[i] = -temp48b[i];
            }

            double* cdea = stackalloc double[1152];
            int cdealen = ArithmeticFunctions.FastExpansionSumZeroelim(temp48alen, temp48a,
                                                  temp48blen, temp48b, cdea);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(cdealen, cdea, pb[0], temp192);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(xlen, temp192, pb[0], det384x);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(cdealen, cdea, pb[1], temp192);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(ylen, temp192, pb[1], det384y);
            zlen = ArithmeticFunctions.ScaleExpansionZeroelim(cdealen, cdea, pb[2], temp192);
            zlen = ArithmeticFunctions.ScaleExpansionZeroelim(zlen, temp192, pb[2], det384z);
            xylen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, det384x, ylen, det384y, detxy);
            double* bdet = stackalloc double[1152];
            int blen = ArithmeticFunctions.FastExpansionSumZeroelim(xylen, detxy, zlen, det384z, bdet);

            temp48alen = ArithmeticFunctions.FastExpansionSumZeroelim(eablen, eab, deblen, deb, temp48a);
            temp48blen = ArithmeticFunctions.FastExpansionSumZeroelim(abdlen, abd, dealen, dea, temp48b);
            for (int i = 0; i < temp48blen; i++)
            {
                temp48b[i] = -temp48b[i];
            }

            double* deab = stackalloc double[96];
            int deablen = ArithmeticFunctions.FastExpansionSumZeroelim(temp48alen, temp48a,
                                                  temp48blen, temp48b, deab);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(deablen, deab, pc[0], temp192);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(xlen, temp192, pc[0], det384x);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(deablen, deab, pc[1], temp192);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(ylen, temp192, pc[1], det384y);
            zlen = ArithmeticFunctions.ScaleExpansionZeroelim(deablen, deab, pc[2], temp192);
            zlen = ArithmeticFunctions.ScaleExpansionZeroelim(zlen, temp192, pc[2], det384z);
            xylen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, det384x, ylen, det384y, detxy);
            double* cdet = stackalloc double[1152];
            int clen = ArithmeticFunctions.FastExpansionSumZeroelim(xylen, detxy, zlen, det384z, cdet);

            temp48alen = ArithmeticFunctions.FastExpansionSumZeroelim(abclen, abc, eaclen, eac, temp48a);
            temp48blen = ArithmeticFunctions.FastExpansionSumZeroelim(bcelen, bce, eablen, eab, temp48b);
            for (int i = 0; i < temp48blen; i++)
            {
                temp48b[i] = -temp48b[i];
            }

            double* eabc = stackalloc double[96];
            int eabclen = ArithmeticFunctions.FastExpansionSumZeroelim(temp48alen, temp48a,
                                                  temp48blen, temp48b, eabc);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(eabclen, eabc, pd[0], temp192);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(xlen, temp192, pd[0], det384x);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(eabclen, eabc, pd[1], temp192);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(ylen, temp192, pd[1], det384y);
            zlen = ArithmeticFunctions.ScaleExpansionZeroelim(eabclen, eabc, pd[2], temp192);
            zlen = ArithmeticFunctions.ScaleExpansionZeroelim(zlen, temp192, pd[2], det384z);
            xylen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, det384x, ylen, det384y, detxy);
            double* ddet = stackalloc double[1152];
            int dlen = ArithmeticFunctions.FastExpansionSumZeroelim(xylen, detxy, zlen, det384z, ddet);

            temp48alen = ArithmeticFunctions.FastExpansionSumZeroelim(bcdlen, bcd, abdlen, abd, temp48a);
            temp48blen = ArithmeticFunctions.FastExpansionSumZeroelim(cdalen, cda, abclen, abc, temp48b);
            for (int i = 0; i < temp48blen; i++)
            {
                temp48b[i] = -temp48b[i];
            }

            double* abcd = stackalloc double[96];
            int abcdlen = ArithmeticFunctions.FastExpansionSumZeroelim(temp48alen, temp48a,
                                                  temp48blen, temp48b, abcd);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(abcdlen, abcd, pe[0], temp192);
            xlen = ArithmeticFunctions.ScaleExpansionZeroelim(xlen, temp192, pe[0], det384x);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(abcdlen, abcd, pe[1], temp192);
            ylen = ArithmeticFunctions.ScaleExpansionZeroelim(ylen, temp192, pe[1], det384y);
            zlen = ArithmeticFunctions.ScaleExpansionZeroelim(abcdlen, abcd, pe[2], temp192);
            zlen = ArithmeticFunctions.ScaleExpansionZeroelim(zlen, temp192, pe[2], det384z);
            xylen = ArithmeticFunctions.FastExpansionSumZeroelim(xlen, det384x, ylen, det384y, detxy);
            double* edet = stackalloc double[1152];
            int elen = ArithmeticFunctions.FastExpansionSumZeroelim(xylen, detxy, zlen, det384z, edet);
            double* abdet = stackalloc double[2304];
            int ablen = ArithmeticFunctions.FastExpansionSumZeroelim(alen, adet, blen, bdet, abdet);
            double* cddet = stackalloc double[2304];
            int cdlen = ArithmeticFunctions.FastExpansionSumZeroelim(clen, cdet, dlen, ddet, cddet);
            double* cdedet = stackalloc double[3456];
            cdelen = ArithmeticFunctions.FastExpansionSumZeroelim(cdlen, cddet, elen, edet, cdedet);
            double* deter = stackalloc double[5760];
            int deterlen = ArithmeticFunctions.FastExpansionSumZeroelim(ablen, abdet, cdelen, cdedet, deter);

            return deter[deterlen - 1];
        }
    }
}
