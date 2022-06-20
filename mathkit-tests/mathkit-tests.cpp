#include "CppUnitTest.h"

import <iostream>;
import mathkit;

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace mathkittests
{
	TEST_CLASS(mathkittests)
	{
	public:
		
        TEST_METHOD(TestMethod1)
        {
            for (int rows = 2; rows < 10; rows++)
            {
                for (int columns = 2; columns < 10; columns++)
                {
                    for (int count = 0; count < 500; count++)
                    {
                        RingMatrix m = RingMatrix::createRandomRingMatrix(rows, rows, 4);
                        std::tuple<RingMatrix, RingMatrix, std::vector<int>> hnf = m.HermiteNormalForm();
                        const RingMatrix U = std::get<0>(hnf);
                        const RingMatrix A = std::get<1>(hnf);
                        if (!(U * m == A))
                        {
                            //std::cout << m;
                            //std::cout << U * m;
                            //std::cout << U;
                            //std::cout << A;
                        }
                        const RingMatrix result = U * m;

                        Assert::IsTrue(U * m == A, L"Failed");
                    }
                }
            }
        }

        

		TEST_METHOD(TestMethod2)
		{
            for (int rows = 2; rows < 10; rows++)
            {
                for (int columns = 2; columns < 10; columns++)
                {
                    for (int count = 0; count < 500; count++)
                    {
                        RingMatrix m = RingMatrix::createRandomRingMatrix(rows, columns, 4);
                        std::tuple<RingMatrix, RingMatrix, RingMatrix> smf = m.SmithNormalForm();
                        const RingMatrix U = std::get<0>(smf);
                        const RingMatrix V = std::get<1>(smf);
                        const RingMatrix A = std::get<2>(smf);
                        if (!(U * m * V == A))
                        {
                            //std::cout << U * m * V;
                            //std::cout << A;
                        }
                        Assert::IsTrue(U * m * V == A, L"Failed");
                    }
                }
            }
           
		}
	};
}
