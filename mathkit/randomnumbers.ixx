export module randomnumbers;

import <random>;

import double3;

export class RandomNumber
{
public:
    static double Uniform()
    {
        return getInstance().uniformDistribution_(getInstance().mt);
    }
    static double Gaussian()
    {
        return getInstance().normalDistribution_(getInstance().mt);
    }
private:
    RandomNumber()
    {
        std::random_device rd;
        //mt = std::mt19937_64(rd());
        mt = std::mt19937_64(14);
        uniformDistribution_ = std::uniform_real_distribution<double>(0.0, 1.0);
        normalDistribution_ = std::normal_distribution<double>();
    }
    ~RandomNumber() {}
    
    static RandomNumber& getInstance() {
        static RandomNumber s;
        return s;
    }

    RandomNumber(RandomNumber const&) = delete;
    RandomNumber& operator= (RandomNumber const&) = delete;

    static inline double3 UnitSphere()
    {
        double ran1, ran2, ranh, ransq;

        do
        {
            ran1 = 2.0 * RandomNumber::Uniform();
            ran2 = 2.0 * RandomNumber::Uniform();
            ransq = ran1 * ran1 + ran2 * ran2;
        } while (ransq >= 1.0);

        ranh = 2.0 * std::sqrt(1.0 - ransq);

        return double3(ran1 * ranh, ran2 * ranh, 1.0 - 2.0 * ransq);
    }

    std::mt19937_64 mt;
    std::uniform_real_distribution<double> uniformDistribution_;
    std::normal_distribution<double> normalDistribution_;
};
