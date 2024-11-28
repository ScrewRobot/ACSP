#include "doctest/doctest.h"
#include <iostream>
#include "ACSP.hpp"
using namespace FastMath;
using namespace ACSP::LTI;
using namespace ACSP::Controller;


TEST_CASE("LADRC2 Test : basic work") {
    LADRC2 adrc;
	adrc.wc = 10;
	adrc.wo = 100;
	CHECK(adrc.b0 == 1.0);
	CHECK(adrc.bound[0] == -100);
	CHECK(adrc.bound[1] == 100);

	constexpr double dt = 1E-3;
	SUBCASE("Check disable")
	{
		adrc.Target = 1;
		adrc.Feedback = 0;
		adrc.Enable = false;
		for (int k = 0; k < 10; ++k)
		{
			adrc.step(dt);
//		cout << adrc.ControlOutput << endl;
		}
		CHECK(adrc.ControlOutput == 0);

	}
	SUBCASE("Check enable")
	{
		adrc.Target = 1;
		adrc.TargetSpeed = 0;
		adrc.Feedback = 0;
		adrc.Enable = true;
		adrc.resetESO(adrc.Feedback);

//		cout << "before step " << endl;
//		cout << adrc.ESO.x << endl;

		for (int k = 0; k < 1000; ++k)
		{
			adrc.step(dt);
//			cout << "step "  << k << endl << adrc.ESO.x(0) << endl;
		}
		CHECK(adrc.ControlOutput == adrc.bound[1]);
	}

	SUBCASE("Check Control direction case 1")
	{
		adrc.Target = -1;
		adrc.TargetSpeed = 0;
		adrc.Feedback = 0;
		adrc.Enable = true;
		adrc.resetESO(adrc.Feedback);

		for (int k = 0; k < 1000; ++k)
		{
			adrc.step(dt);
		}
		CHECK(adrc.ControlOutput == adrc.bound[0]);
	}

	SUBCASE("Check Control direction case 2")
	{
		adrc.Target = 0;
		adrc.TargetSpeed = 0;
		adrc.Feedback = 1;
		adrc.Enable = true;
		adrc.resetESO(adrc.Feedback);

		for (int k = 0; k < 1000; ++k)
		{
			adrc.step(dt);
		}
		CHECK(adrc.ControlOutput == adrc.bound[0]);
	}

	SUBCASE("Check Control direction case 3")
	{
		adrc.Target = 0;
		adrc.TargetSpeed = 0;
		adrc.Feedback = -10;
		adrc.Enable = true;
		adrc.resetESO(adrc.Feedback);

		for (int k = 0; k < 1000; ++k)
		{
			adrc.step(dt);
		}
		CHECK(adrc.ControlOutput == adrc.bound[1]);
	}
}

TEST_CASE("LADRC2 Test : basic control")
{

	using namespace std;

    StateSpace<2, 1, 1> Plant;
    Plant.cleanAll();
    Plant.A(0, 1) = 1;
    Plant.B(1, 0) = 1;
    Plant.C(0, 0) = 1;

    LADRC2 adrc;
    adrc.setParam(10, 100, 1, {-100, 100});

	CHECK(adrc.wc == 10);
	CHECK(adrc.wo == 100);
	CHECK(adrc.b0 == 1.0);
	CHECK(adrc.bound[0] == -100);
	CHECK(adrc.bound[1] == 100);

	adrc.resetESO(Plant.y(0));
	adrc.Enable = false;

	SUBCASE("Target 1")
	{
		double log_maxoutput = 0;
		double log_minoutput = 0;
		for (int k = 0; k < 10000; ++k)
		{
			if (k > 10)
				adrc.Enable = true;

			adrc.Target = 1;
			adrc.TargetSpeed = 0;
			adrc.Feedback = Plant.y(0);
			Plant.u(0) = adrc.ControlOutput;

			adrc.step(1E-3);
			Plant.step(1E-3);

			log_maxoutput = (log_maxoutput<adrc.ControlOutput)?adrc.ControlOutput:log_maxoutput;
			log_minoutput = (log_minoutput>adrc.ControlOutput)?adrc.ControlOutput:log_minoutput;

//		cout << "t : " << Plant.t << " Plant : " << Plant.y(0) << " adrc out: " << adrc.ControlOutput << endl;
		}
		CHECK(Plant.y(0) == doctest::Approx(1).epsilon(1E-6));
		CHECK(adrc.ControlOutput == doctest::Approx(0).epsilon(1E-6));
		CHECK(log_maxoutput <= adrc.bound[1]);
		CHECK(log_minoutput >= adrc.bound[0]);
	}

	SUBCASE("Target -10")
	{
		double log_maxoutput = 0;
		double log_minoutput = 0;
		for (int k = 0; k < 10000; ++k)
		{
			if (k > 10)
				adrc.Enable = true;

			adrc.Target = -10;
			adrc.TargetSpeed = 0;
			adrc.Feedback = Plant.y(0);
			Plant.u(0) = adrc.ControlOutput;

			adrc.step(1E-3);
			Plant.step(1E-3);

			log_maxoutput = (log_maxoutput<adrc.ControlOutput)?adrc.ControlOutput:log_maxoutput;
			log_minoutput = (log_minoutput>adrc.ControlOutput)?adrc.ControlOutput:log_minoutput;

//		cout << "t : " << Plant.t << " Plant : " << Plant.y(0) << " adrc out: " << adrc.ControlOutput << endl;
		}
		CHECK(Plant.y(0) == doctest::Approx(-10).epsilon(1E-6));
		CHECK(adrc.ControlOutput == doctest::Approx(0).epsilon(1E-6));
		CHECK(log_maxoutput <= adrc.bound[1]);
		CHECK(log_minoutput >= adrc.bound[0]);
	}

}

TEST_CASE("LADRC2 Test : spring dagger")
{
	using namespace LTI;
	using namespace std;

	double k = 1;
	double m = 0.1;
	double B = 0.2;



	StateSpace<2, 1, 1> Plant;

    Plant.cleanAll();
    Plant.A(0, 1) = 1;
    Plant.A(1, 0) = -k/m;
    Plant.A(1, 1) = -B/m;
    Plant.B(1, 0) = 1/m;
    Plant.C(0, 0) = 1;

	LADRC2 adrc;
	// LADRC-DesignTools wf=20 rad/s, PM=60 deg
	adrc.setParam(5.0326, 79.4823, 17.0905, {-100, 100});

	adrc.resetESO(Plant.y(0));
	adrc.Enable = false;

	SUBCASE("target = 1")
	{
		double log_maxoutput = 0;
		double log_minoutput = 0;
		for (int k = 0; k < 10000; ++k)
		{
			if (k > 10)
				adrc.Enable = true;

			adrc.Target = 1;
			adrc.TargetSpeed = 0;
			adrc.Feedback = Plant.y(0);
			Plant.u(0) = adrc.ControlOutput;

			adrc.step(1E-3);
			Plant.step(1E-3);

			log_maxoutput = (log_maxoutput<adrc.ControlOutput)?adrc.ControlOutput:log_maxoutput;
			log_minoutput = (log_minoutput>adrc.ControlOutput)?adrc.ControlOutput:log_minoutput;

//	cout << "t : " << Plant.t << ", Plant : " << Plant.y(0) << ", adrc out: " << adrc.ControlOutput << endl;
		}
		CHECK(Plant.y(0) == doctest::Approx(1).epsilon(1E-6));
		CHECK(adrc.ControlOutput == doctest::Approx(1).epsilon(1E-6));
		CHECK(log_maxoutput <= adrc.bound[1]);
		CHECK(log_minoutput >= adrc.bound[0]);
	}

	SUBCASE("target = -10")
	{
		double log_maxoutput = 0;
		double log_minoutput = 0;
		for (int k = 0; k < 10000; ++k)
		{
			if (k > 10)
				adrc.Enable = true;

			adrc.Target = -10;
			adrc.TargetSpeed = 0;
			adrc.Feedback = Plant.y(0);
			Plant.u(0) = adrc.ControlOutput;

			adrc.step(1E-3);
			Plant.step(1E-3);

			log_maxoutput = (log_maxoutput<adrc.ControlOutput)?adrc.ControlOutput:log_maxoutput;
			log_minoutput = (log_minoutput>adrc.ControlOutput)?adrc.ControlOutput:log_minoutput;

//	cout << "t : " << Plant.t << ", Plant : " << Plant.y(0) << ", adrc out: " << adrc.ControlOutput << endl;
		}
		CHECK(Plant.y(0) == doctest::Approx(-10).epsilon(1E-6));
	}

}


TEST_CASE("LADRC2 Test : reset during run")
{
	using namespace LTI;
	using namespace std;

	double k = 1;
	double m = 0.1;
	double B = 0.2;

	StateSpace<2, 1, 1> Plant;
    Plant.cleanAll();
    Plant.A(0, 1) = 1;
    Plant.A(1, 0) = -k/m;
    Plant.A(1, 1) = -B/m;
    Plant.B(1, 0) = 1/m;
    Plant.C(0, 0) = 1;

	LADRC2 adrc;
	// LADRC-DesignTools wf=20 rad/s, PM=60 deg
	adrc.setParam(5.0326, 79.4823, 17.0905, {-100, 100});

	adrc.resetESO(Plant.y(0));
	adrc.Enable = false;


	double log_maxoutput = 0;
	double log_minoutput = 0;
	for (int k = 0; k < 10000; ++k)
	{
		if (k > 10)
			adrc.Enable = true;

		if (Plant.t == 2)
			adrc.resetESO(Plant.y(0));
		if (Plant.t == 3)
			adrc.resetESO(Plant.y(0));
		adrc.Target = 1;
		adrc.TargetSpeed = 0;
		adrc.Feedback = Plant.y(0);
		Plant.u(0) = adrc.ControlOutput;

		adrc.step(1E-3);
		Plant.step(1E-3);

		log_maxoutput = (log_maxoutput<adrc.ControlOutput)?adrc.ControlOutput:log_maxoutput;
		log_minoutput = (log_minoutput>adrc.ControlOutput)?adrc.ControlOutput:log_minoutput;

//	cout << "t : " << Plant.t << ", Plant : " << Plant.y(0) << ", adrc out: " << adrc.ControlOutput << ", ESO : " << adrc.ESO.x(0) << endl;
	}
	CHECK(Plant.y(0) == doctest::Approx(1).epsilon(1E-6));
	CHECK(adrc.ControlOutput == doctest::Approx(1).epsilon(1E-6));
	CHECK(log_maxoutput <= adrc.bound[1]);
	CHECK(log_minoutput >= adrc.bound[0]);


}

TEST_CASE(" ControllabilityMatrix and  test")
{
	using namespace FastMath;
    SquareMatrix<double, 3> A({0, 1, 3,
                               -10, -2, 3,
                               1, 2, 3});
    Matrix<double, 3, 2> B({0, 1, 10, 2, 5, 3});
    Matrix<double, 1, 3> C;
    C(0,0) = 1;
    StateSpace<3, 2, 1> Plant;
    Plant.A = A;
    Plant.B = B;
    Plant.C = C;

    auto Pc = Plant.ControllabilityMatrix();
    Matrix<double, 3, 6> Pc_matlab({0  ,  1  ,  25  ,  11  , 100  ,  37,
                                    10 ,   2 ,   -5 ,   -5 , -135 ,  -58,
                                    5  ,  3  ,  35  ,  14  , 120  ,  43});
    CHECK(Pc == Pc_matlab);


    auto Po = Plant.ObservabilityMatrix();
    Matrix<double, 3, 3> Po_matlab({     1 ,    0 ,    0,
                                        0  ,   1  ,   3,
                                        -7 ,    4 ,   12});
    CHECK(Po == Po_matlab);
}


TEST_CASE(" C2D and  test")
{
    SquareMatrix<double, 3> A({1,2,3,4,5,6,7,8,9});
    Matrix<double, 3, 1> B({0,0,1});
    Matrix<double, 1, 3> C({1, 0, 0});

    StateSpace css(A,B,C);

    auto dss1 = C2D(css, 1E-3, C2D_TYPE::Tustin);
    auto dss2 = C2D(css, 1E-3, C2D_TYPE::ZOH);

    Vector<double, 3> x({1.0,2.0,3.0});
    Vector<double, 1> u({10.0});

    css.setState(x);
    dss1.setState(x);
    dss2.setState(x);

    css.setInput(u);
    dss1.setInput(u);
    dss2.setInput(u);

    css.step(1e-3);
    dss1.step(1e-3);
    dss2.step(1e-3);

    auto res = css.getOutput();
    auto res1 = dss1.getOutput();
    auto res2 = dss2.getOutput();

    double tustin_res = 1.0213261180077715106762070718105;  //matlab tustin result
    double zoh_res = 1.0141296847598906261112006177427;     //matlab zoh result

    CHECK(res1(0) == doctest::Approx(tustin_res).epsilon(1E-6));
    CHECK(res2(0) == doctest::Approx(zoh_res).epsilon(1E-6));


}
