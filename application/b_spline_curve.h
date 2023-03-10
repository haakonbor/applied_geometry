#ifndef B_SPLINE_CURVE_H
#define B_SPLINE_CURVE_H

#include <numeric>
#include <parametrics/gmpcurve.h>

namespace hbb {

template <typename T> class My_B_spline : public GMlib::PCurve<T,3> {
    GM_SCENEOBJECT(My_B_spline)

public:
    My_B_spline(const GMlib::DVector<GMlib::Vector<float,3>>& c);
    My_B_spline(const GMlib::DVector<GMlib::Vector<float,3>>& p, int n);
    virtual ~My_B_spline();

    //****************************************
    //****** Virtual public functions   ******
    //****************************************

    // from PCurve
    bool isClosed() const override;

protected:
    // Virtual functions from PCurve, which have to be implemented locally
    void eval(T t, int d, bool l) const override;
    T getStartP() const override;
    T getEndP() const override;

    void createKnotVector(T start, T end);
    T getW(int d, int i, T knot_val) const;
    int getKnotIndex(T t) const;
    GMlib::Vector<T,3> getBasisVector(int i, T t) const;

    GMlib::DVector<GMlib::Vector<T, 3>> _control_points;
    std::vector<T> _knot_vector;

    int _k = 3;     // Order
    int _n;         // Number of control points
    int _d = 2;     // Degree



}; // END class My_B_spline

//*****************************************
// Constructors and destructor           **
//*****************************************

template <typename T>
inline My_B_spline<T>::My_B_spline(const GMlib::DVector<GMlib::Vector<float,3>>& c)
    : GMlib::PCurve<T, 3>(0,0,0), _control_points(c), _n(c.getDim()) {
    createKnotVector(T(0), T(1));
    std::cout << "Control points: " << _control_points << std::endl;
}

template <typename T>
inline My_B_spline<T>::My_B_spline(const GMlib::DVector<GMlib::Vector<float,3>>& p, int n)
    : GMlib::PCurve<T, 3>(0,0,0), _n(n) {
    // From formula 6.21, 6.22, and 6.23 in chapter 6.5 in the book
    auto x = std::vector<T>{T(0)};
    for (int i = 1; i < p.getDim(); i++) {
        x.push_back((p[i] - p[i-1]).getLength() + x.back());
    }

    createKnotVector(x[0], x.back());
    std::cout << "X: " << x << std::endl;

    auto A = GMlib::DMatrix<T>(p.getDim(), n, T(0));
    for (int i = 0; i < x.size(); i++) {
        int knot_index = getKnotIndex(x[i]);
        auto basis_vector = getBasisVector(knot_index, x[i]);
        A[i][knot_index - 2] = basis_vector[0];
        A[i][knot_index - 1] = basis_vector[1];
        A[i][knot_index]     = basis_vector[2];
    }

    auto AT = A.transpose();
    auto B = AT * A;
    auto y = AT * p;
    B.invert();
    _control_points = B * y;

    std::cout << "Control points: " << _control_points << std::endl;
}

template <typename T>
My_B_spline<T>::~My_B_spline() {}

//**************************************
//        Public local functons       **
//**************************************



//***************************************************
// Overrided (public) virtual functons from PCurve **
//***************************************************

template <typename T>
bool My_B_spline<T>::isClosed() const { return false; }

template <typename T> void My_B_spline<T>::eval(T t, int d, bool /*l*/) const {

    this->_p.setDim(d + 1);

    int i = getKnotIndex(t);

    auto basis = getBasisVector(i, t);

    this->_p[0] = basis[0] * _control_points[i-2] + basis[1] * _control_points[i-1] + basis[2] * _control_points[i];
}

template <typename T> T My_B_spline<T>::getStartP() const { return _knot_vector[_k]; }

template <typename T> T My_B_spline<T>::getEndP() const { return _knot_vector[_knot_vector.size() - _k]; }


//**************************************
//        Protected local functons    **
//**************************************

template <typename T>
void My_B_spline<T>::createKnotVector(T start, T end) {
    _knot_vector.clear();

    // Since we want a clamped B-spline curve, we set the first k and last k knot values to be equal
    for (int i = 0; i < _k; i++) {
        _knot_vector.push_back(start);
    }

    for (int i = 1; i <= _n -_k; i++) {
        _knot_vector.push_back(start + i * (end - start) / (_n - _k + 1));
    }

    for (int i = 0; i < _k; i++) {
        _knot_vector.push_back(end);
    }

    // Debug
    std::cout << "NUMBER OF CONTROL POINTS: " << _n << std::endl;
    std::cout << "KNOT VECTOR: " << _knot_vector << std::endl;
}

template <typename T>
T My_B_spline<T>::getW(int d, int i, T t) const {
    // From function 6.11 at page 82 in the book
    return (t - _knot_vector[i]) / (_knot_vector[i+d] - _knot_vector[i]);
}

template <typename T>
GMlib::Vector<T,3> My_B_spline<T>::getBasisVector(int i, T t) const {
    /*
     * c(t) = ( 1 - w_1_i(t)   w_1_i(t) ) ( 1 - w_2_i-1(t)   w_2_i-1(t)          0     ) ( c_i-2 )
     *                                    (       0         1 - w_2_i(t)      w_2_i(t) ) ( c_i-1 )
     *                                                                                   (  c_i  )
    */

    auto w_1_i = getW(1, i, t);
    auto w_2_i = getW(2, i, t);
    auto w_2_i_1 = getW(2, i - 1, t);

    auto basis1 = (1 - w_1_i) * (1 - w_2_i_1);
    auto basis2 = (1 - w_1_i) * w_2_i_1 + w_1_i * (1 - w_2_i);
    auto basis3 = w_1_i * w_2_i;

    return {basis1, basis2, basis3};
}

template <typename T>
int My_B_spline<T>::getKnotIndex(T t) const {
    // Since the curve is open/clamped: t_d <= t <= t_n
    for (int i = _d; i < _n; i++) {
        if (t <= _knot_vector[i + 1]){
            return i;
        }
    }
}

}; // END namespace hbb

#endif // B_SPLINE_CURVE_H
