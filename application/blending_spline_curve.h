#ifndef BLENDING_SPLINE_CURVE_H
#define BLENDING_SPLINE_CURVE_H

#include <numeric>
#include <parametrics/gmpcurve.h>
#include <parametrics/curves/gmpsubcurve.h>

namespace hbb {

template <typename T> class Blending_spline_curve : public GMlib::PCurve<T,3> {
    GM_SCENEOBJECT(Blending_spline_curve)

public:
    Blending_spline_curve(GMlib::PCurve<T,3>* model_curve, int n);
    virtual ~Blending_spline_curve();

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
    void createLocalCurves();
    T getW(int d, int i, T knot_val) const;
    int getKnotIndex(T t) const;
    GMlib::Vector<T,3> getBasisVector(int i, T t) const;
    float bFunction(T t) const;

    void localSimulate(double dt) override;

    GMlib::DVector<GMlib::Vector<T, 3>> _control_points;
    std::vector<T> _knot_vector;
    GMlib::PCurve<T, 3>* _model_curve;
    std::vector<GMlib::PSubCurve<T>*> _local_curves;

    int _n;         // Number of control points
    int _k = 2;

    int _count = 0;
    int _direction = 1;
    float _x = 0.0f;

    bool _is_closed;



}; // END class Blending_spline_curve

//*****************************************
// Constructors and destructor           **
//*****************************************

template <typename T>
inline Blending_spline_curve<T>::Blending_spline_curve(GMlib::PCurve<T, 3>* model_curve, int n)
    : GMlib::PCurve<T, 3>(0,0,0), _model_curve(model_curve), _n(n), _is_closed(model_curve->isClosed()) {

    createKnotVector(model_curve->getParStart(), model_curve->getParEnd());
    createLocalCurves();
}

template <typename T>
Blending_spline_curve<T>::~Blending_spline_curve() {}

//**************************************
//        Public local functons       **
//**************************************



//***************************************************
// Overrided (public) virtual functons from PCurve **
//***************************************************

template <typename T>
bool Blending_spline_curve<T>::isClosed() const { return _is_closed; }

template <typename T> void Blending_spline_curve<T>::eval(T t, int d, bool /*l*/) const {

    this->_p.setDim(d + 1);

    int i = getKnotIndex(t);

    auto basis = getBasisVector(i, t);

    this->_p = basis[0] * _local_curves[i-1]->evaluateParent(t, d) + basis[1] * _local_curves[i]->evaluateParent(t, d);
}

template <typename T> T Blending_spline_curve<T>::getStartP() const { return _model_curve->getParStart(); }

template <typename T> T Blending_spline_curve<T>::getEndP() const { return _model_curve->getParEnd(); }


//**************************************
//        Protected local functons    **
//**************************************

template <typename T>
void Blending_spline_curve<T>::createKnotVector(T start, T end) {
    _knot_vector.clear();

    _knot_vector.push_back(start);

    for (int i = 0; i < _n + 1; i++) {
        _knot_vector.push_back(start + i * (_model_curve->getParDelta() / _n));
    }

    _knot_vector[0] = _knot_vector[1] - _knot_vector[_n + 1] + _knot_vector[_n];    // Knot vetor bites itself in the tail

    // Debug
    std::cout << "KNOT VECTOR: " << _knot_vector << std::endl;
}

template <typename T>
void Blending_spline_curve<T>::createLocalCurves() {
    for (int i = 0; i < _n; i++) {
        // subcurve(original_curve, start_parameter, end_parameter, parameter_local_origin)
        _local_curves.push_back(new GMlib::PSubCurve<T>(_model_curve, _knot_vector[i], _knot_vector[i+2], _knot_vector[i+1]));
        _local_curves[i]->toggleDefaultVisualizer();
        _local_curves[i]->setCollapsed(true);
        _local_curves[i]->setParent(this);
    }

    _local_curves.push_back(_local_curves[0]);
}

template <typename T>
T Blending_spline_curve<T>::getW(int d, int i, T t) const {
    // From function 6.11 at page 82 in the book
    return (t - _knot_vector[i]) / (_knot_vector[i+d] - _knot_vector[i]);
}

template <typename T>
GMlib::Vector<T,3> Blending_spline_curve<T>::getBasisVector(int i, T t) const {
    T w_1_i = getW(1, i, t);
    T B = bFunction(w_1_i);

    return GMlib::Vector<T, 2> {1-B, B};
}

template <typename T>
float Blending_spline_curve<T>::bFunction(T t) const {
    // Symmetric TB-function of order 2 (page 134 in the book)
    return t - T(1 / M_2PI) * sin(M_2PI * t);
}

template <typename T>
int Blending_spline_curve<T>::getKnotIndex(T t) const {
    if (t >= getEndP()) {
        return _knot_vector.size() - 2;
    }

    return std::distance(_knot_vector.begin(), std::upper_bound(_knot_vector.begin(), _knot_vector.end(), t)) - 1;
}

template <typename T>
void Blending_spline_curve<T>::localSimulate(double dt) {
    for (int i = 0; i < _n; i++) {
        if (i % 2 == 0) {
            _local_curves[i]->translate({T(0.1) * cos(_x), T(0.1) * sin(_x), T(0.0)});
        }
        else {
            _local_curves[i]->rotate(cos(_x) * 0.1, {1.0, 0.0, 0.0});
        }
    }

    _x += 0.1f;
    if (_x >= M_2PI) {
        _x = 0.0f;
    }


    this->resample();
    this->setEditDone();

}

}; // END namespace hbb

#endif // BLENDING_SPLINE_CURVE_H
