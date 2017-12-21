using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace CoreMath
{
    public enum EulerOrder
    {
        XYZ,
        YXZ,
        ZXY,
        ZYX,
        YZX,
        XZY
    }

    public static class EulerExtensions
    {
        public static float[] EulerFromQuaternion(this float[] quaternion, EulerOrder order)
        {
            var matrix = quaternion.MatrixFromQuaternion();

            return matrix.EulerFromMatrix(order);
        }

        public static float[] EulerFromMatrix(this float[] matrix, EulerOrder order)
        {
            var result = new float[] { 0, 0, 0 };
            // assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)
            var m11 = matrix[0]; var m12 = matrix[4]; var m13 = matrix[8];
            var m21 = matrix[1]; var m22 = matrix[5]; var m23 = matrix[9];
            var m31 = matrix[2]; var m32 = matrix[6]; var m33 = matrix[10];

            if (order == EulerOrder.XYZ)
            {

                result[1] = (float)Math.Asin(m13.Clamp(-1, 1));

                if (Math.Abs(m13) < 0.99999)
                {

                    result[0] = (float)Math.Atan2(-m23, m33);
                    result[2] = (float)Math.Atan2(-m12, m11);

                }
                else
                {

                    result[0] = (float)Math.Atan2(m32, m22);
                    result[2] = 0;

                }

            }
            else if (order == EulerOrder.YXZ)
            {

                result[0] = (float)Math.Asin(-(m23.Clamp(-1, 1)));

                if (Math.Abs(m23) < 0.99999)
                {

                    result[1] = (float)Math.Atan2(m13, m33);
                    result[2] = (float)Math.Atan2(m21, m22);

                }
                else
                {

                    result[1] = (float)Math.Atan2(-m31, m11);
                    result[2] = 0;

                }

            }
            else if (order == EulerOrder.ZXY)
            {

                result[0] = (float)Math.Asin(m32.Clamp(-1, 1));

                if (Math.Abs(m32) < 0.99999)
                {

                    result[1] = (float)Math.Atan2(-m31, m33);
                    result[2] = (float)Math.Atan2(-m12, m22);

                }
                else
                {

                    result[1] = 0;
                    result[2] = (float)Math.Atan2(m21, m11);

                }

            }
            else if (order == EulerOrder.ZYX)
            {

                result[1] = (float)Math.Asin(-(m31.Clamp(-1, 1)));

                if (Math.Abs(m31) < 0.99999)
                {

                    result[0] = (float)Math.Atan2(m32, m33);
                    result[2] = (float)Math.Atan2(m21, m11);

                }
                else
                {

                    result[0] = 0;
                    result[2] = (float)Math.Atan2(-m12, m22);

                }

            }
            else if (order == EulerOrder.YZX)
            {

                result[2] = (float)Math.Asin(m21.Clamp(-1, 1));

                if (Math.Abs(m21) < 0.99999)
                {

                    result[0] = (float)Math.Atan2(-m23, m22);
                    result[1] = (float)Math.Atan2(-m31, m11);

                }
                else
                {

                    result[0] = 0;
                    result[1] = (float)Math.Atan2(m13, m33);

                }

            }
            else if (order == EulerOrder.XZY)
            {

                result[2] = (float)Math.Asin(-(m12.Clamp(-1, 1)));

                if (Math.Abs(m12) < 0.99999)
                {

                    result[0] = (float)Math.Atan2(m32, m22);
                    result[1] = (float)Math.Atan2(m13, m11);

                }
                else
                {

                    result[0] = (float)Math.Atan2(-m23, m33);
                    result[1] = 0;

                }

            }

            return result;
        }


    }
}
