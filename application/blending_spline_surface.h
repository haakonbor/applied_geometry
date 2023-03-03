#ifndef BLENDING_SPLINE_SURFACE_H
#define BLENDING_SPLINE_SURFACE_H

#include "simplesubsurf.h"

namespace hbb {

template <typename T> class Blending_spline_surface : public GMlib::PSurf<T, 3> {
    GM_SCENEOBJECT(Blending_spline_surface)

public:
    Blending_spline_surface(GMlib::PSurf<T, 3>* model_surface, int n_u = 4, int n_v = 4, int d = 1);
    virtual ~Blending_spline_surface();

    bool isClosedU() const override;
    bool isClosedV() const override;

protected:
    T    getStartPU() const override;
    T    getEndPU() const override;
    T    getStartPV() const override;
    T    getEndPV() const override;
    void eval(T u, T v, int d1, int d2, bool lu = true,
              bool lv = true) const override;
    void localSimulate(double dt) override;

private:
    int                                 _n_u, _n_v, _k, _d;
    bool                                _is_U_closed, _is_V_closed;
    GMlib::PSurf<T, 3>*                 _model_surface;
    std::vector<T>                      _u;
    std::vector<T>                      _v;
    GMlib::DMatrix<GMlib::PSurf<T, 3>*> _s;

    void           createU(T start, T end);
    void           createV(T start, T end);
    void           createLocalSurfaces();
    int            getIndex(T t, const std::vector<T>& _t, bool is_closed, int n) const;
    void           getBasis(T& B, T& d_B, int i, T t, const std::vector<T>& _t) const;
    std::vector<T> bFunction(T t) const;
    std::vector<T> _w(int d, int i, T t, const std::vector<T>& _t) const;

}; // END class Blending_spline_surface

template <typename T>
inline Blending_spline_surface<T>::Blending_spline_surface(GMlib::PSurf<T, 3>* model_surface,
                                                           int n_u, int n_v, int d)
    : _model_surface(model_surface), _n_u(n_u), _n_v(n_v), _k(d + 1), _d(d),
    _is_U_closed(model_surface->isClosedU()), _is_V_closed(model_surface->isClosedV()) {

    _s.setDim((isClosedU() ? _n_u + 1 : _n_u),
              (isClosedV() ? _n_v + 1 : _n_v));

    createU(model_surface->getParStartU(), model_surface->getParEndU());
    createV(model_surface->getParStartV(), model_surface->getParEndV());

    createLocalSurfaces();
}

template <typename T>
Blending_spline_surface<T>::~Blending_spline_surface() { }

template <typename T>
bool Blending_spline_surface<T>::isClosedU() const {
    return _is_U_closed;
}

template <typename T>
bool Blending_spline_surface<T>::isClosedV() const {
    return _is_V_closed;
}

template <typename T>
void Blending_spline_surface<T>::eval(T u, T v, int d1, int d2,
                                      bool /*lu = true*/,
                                      bool /*lv = true*/) const {
    this->_p.setDim(d1 + 1, d2 + 1);

    int u_i = getIndex(u, _u, isClosedU(), _n_u);
    int v_i = getIndex(v, _v, isClosedV(), _n_v);

    T B_u, d_B_u;
    T B_v, d_B_v;
    getBasis(B_u, d_B_u, u_i, u, _u);
    getBasis(B_v, d_B_v, v_i, v, _v);

    // evaluateParent returns |S,     s_u|
    //                        |s_v,  s_uv|
    auto s_0_0 = _s(u_i - 1)(v_i - 1)->evaluateParent(u, v, d1, d2);
    auto s_0_1 = _s(u_i - 1)(v_i)->evaluateParent(u, v, d1, d2);
    auto s_1_0 = _s(u_i)(v_i - 1)->evaluateParent(u, v, d1, d2);
    auto s_1_1 = _s(u_i)(v_i)->evaluateParent(u, v, d1, d2);

    this->_p[0][0] = s_0_0(0)(0) * (1 - B_u) * (1 - B_v)
                     + s_1_0(0)(0) * B_u * (1 - B_v)
                     + s_0_1(0)(0) * (1 - B_u) * B_v
                     + s_1_1(0)(0) * B_u * B_v;

    // PARTIAL DERIVATIVES
    // Chain rule: f'(g(x)) = f'(x) * g(x) + f(x) * g'(x)

    // 1-B, B   |  d(1-B) = d(-B) = -d(B), d(B)

    // u
    this->_p[0][1] = s_0_0(0)(0) * -d_B_u * (1 - B_v) + s_0_0(0)(1) * (1 - B_u) * (1 - B_v)
                     + s_1_0(0)(0) * d_B_u * (1 - B_v) + s_1_0(0)(1) * B_u * (1 - B_v)
                     + s_0_1(0)(0) * -d_B_u * B_v + s_0_1(0)(1) * (1 - B_u) * B_v
                     + s_1_1(0)(0) * d_B_u * B_v + s_1_1(0)(1) * B_u * B_v;

    // v
    this->_p[1][0] = s_0_0(0)(0) * (1 - B_u) * -d_B_v + s_0_0(1)(0) * (1 - B_u) * (1 - B_v)
                     + s_1_0(0)(0) * B_u * -d_B_v + s_1_0(1)(0) * B_u * (1 - B_v)
                     + s_0_1(0)(0) * (1 - B_u) * d_B_v + s_0_1(1)(0) * (1 - B_u) * B_v
                     + s_1_1(0)(0) * B_u * d_B_v + s_1_1(1)(0) * B_u * B_v;
}

template <typename T>
T Blending_spline_surface<T>::getStartPU() const {
    return _model_surface->getParStartU();
}

template <typename T>
T Blending_spline_surface<T>::getEndPU() const {
    return _model_surface->getParEndU();
}

template <typename T>
T Blending_spline_surface<T>::getStartPV() const {
    return _model_surface->getParStartV();
}

template <typename T>
T Blending_spline_surface<T>::getEndPV() const {
    return _model_surface->getParEndV();
}

template <typename T>
void Blending_spline_surface<T>::getBasis(T& B, T& d_B, int i, T t, const std::vector<T>& _t) const {
    auto w_1_i = _w(1, i, t, _t);
    auto b = bFunction(w_1_i[0]);
    B = b[0];
    d_B = b[1] * w_1_i[1];
}

template <typename T>
inline std::vector<T> Blending_spline_surface<T>::bFunction(T t) const {
    // Symmetric TB-function of order 2 (page 134 in the book)
    // + derivative
    return {T(t - T(1 / M_2PI) * sin(M_2PI * t)), T(1 - cos(T(M_2PI) * t))};
}

template <typename T>
inline std::vector<T>
Blending_spline_surface<T>::_w(int d, int i, T t, const std::vector<T>& _t) const {
    // From function 6.11 at page 82 in the book
    // + derivative
    return {(t - _t[i]) / (_t[i + d] - _t[i]), 1 / (_t[i + d] - _t[i])};
}

template <typename T>
inline void Blending_spline_surface<T>::createU(T start, T end) {
    _u.clear();
    // Closed knot vector
    if (isClosedU()) {
        _u.push_back(start);

        for (int i = 0; i < _n_u + _d; i++) {
            _u.push_back(start + i * (_model_surface->getParDeltaU() / _n_u));
        }

        _u[0] = _u[1] - _u[_n_u + 1] + _u[_n_u];
    }
    // Clamped knot vector
    else {
        for (int i = 0; i < _k; i++) {
            _u.push_back(start);
        }

        for (int i = 1; i <= _n_u -_k; i++) {
            _u.push_back(start + i * (end - start) / (_n_u - _k + 1));
        }

        for (int i = 0; i < _k; i++) {
            _u.push_back(end);
        }
    }
}

template <typename T>
inline void Blending_spline_surface<T>::createV(T start, T end) {
    _v.clear();
    // Closed knot vector
    if (isClosedV()) {
        _v.push_back(start);

        for (int i = 0; i < _n_v + _d; i++) {
            _v.push_back(start + i * (_model_surface->getParDeltaV() / _n_v));
        }

        _v[0] = _v[1] - _v[_n_v + 1] + _v[_n_v];
    }
    // Clamped knot vector
    else {
        for (int i = 0; i < _k; i++) {
            _v.push_back(start);
        }

        for (int i = 1; i <= _n_v -_k; i++) {
            _v.push_back(start + i * (end - start) / (_n_v - _k + 1));
        }

        for (int i = 0; i < _k; i++) {
            _v.push_back(end);
        }
    }
}

template <typename T>
inline void Blending_spline_surface<T>::createLocalSurfaces() {
    for (int i = 0; i < _n_u; i++) {
        for (int j = 0; j < _n_v; j++) {
            _s[i][j]
                = new PSimpleSubSurf<T>(_model_surface, _u[i], _u[i + 2], _u[i + 1],
                                        _v[j], _v[j + 2], _v[j + 1]);
            _s[i][j]->toggleDefaultVisualizer();
            _s[i][j]->sample(5, 5, 1, 1);
            _s[i][j]->setCollapsed(true);
            _s[i][j]->setParent(this);
            this->insert(_s[i][j]);
        }
    }

    if (isClosedV()) {
        for (int i = 0; i < _n_u; i++) {
            _s[i][_n_v] = _s[i][0];
        }
    }

    if (isClosedU()) {
        for (int j = 0; j < _n_v; j++) {
            _s[_n_u][j] = _s[0][j];
        }
    }

    if (isClosedV() && isClosedU()){
        _s[_n_u][_n_v] = _s[0][0];
    }
}


template <typename T>
void Blending_spline_surface<T>::localSimulate(double dt) {
    this->sample(_visu[0][0][0], _visu[0][0][1], 1, 1);
}

template <typename T>
int Blending_spline_surface<T>::getIndex(T t, const std::vector<T>& _t,
                                         bool is_closed, int n) const {

    if (is_closed) {
        // t_d <= t < t_n+d
        for (int i = _d; i < n + _d - 1; i++) {
            if (t <= _t[i + 1]) {
                return i;
            }
        }
    }
    else {
        // t_d <= t <= t_n
        for (int i = _d; i < n; i++) {
            if (t <= _t[i + 1]) {
                return i;
            }
        }
    }
}

}   // namespace hbb


#endif // BLENDING_SPLINE_SURFACE_H
