using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Text;

namespace CoreMath.Test
{
    [TestClass]
    public class VectorTests
    {
        [TestMethod]
        public void TestAxis()
        {
            var axis = new float[] { 1, 0, 0 }.VectorAxis();

            var dot1 = axis[0].VectorDotProduct(axis[1]);
            var dot2 = axis[0].VectorDotProduct(axis[2]);
            var dot3 = axis[1].VectorDotProduct(axis[2]);

            Assert.AreEqual(0, dot1);
            Assert.AreEqual(0, dot2);
            Assert.AreEqual(0, dot3);
        }
    }
}
