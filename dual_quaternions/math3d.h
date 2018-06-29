#ifndef CVLAB_MATH3D_H
#define CVLAB_MATH3D_H

#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

typedef unsigned char uint8_t;
typedef unsigned int uint32_t;

namespace cvlab {

   static const double pi = M_PI;
   static const double pi2 = 2.0 * M_PI;
   static const double rad_on_deg = pi / 180.;
   static const double deg_on_rad = 180. / pi;

   struct color_rgb24 {
      uint8_t r, g, b;
      explicit color_rgb24() : r(255), g(255), b(255) {}
      color_rgb24(uint8_t R, uint8_t G, uint8_t B) : r(R), g(G), b(B) {}
   };

   const color_rgb24 white(255,255,255);
   const color_rgb24 black(0,0,0);
   const color_rgb24 grey(128,128,128);
   const color_rgb24 red(255,0,0);
   const color_rgb24 green(0,255,0);
   const color_rgb24 blue(0,0,255);
   const color_rgb24 magenta(255,0,255);
   const color_rgb24 cyan(0,255,255);
   const color_rgb24 yellow(255,255,0);
   const color_rgb24 orange(255,128,0);

   color_rgb24 rand_color();

   template <typename T>
   inline bool almost_zero(T a, double e)
   {
      return (a == T(0)) || (a > 0 && a < e) || (a < 0 && a > -e);
   }


   template<typename T>
   struct vec2d
   {
      T x,y;

      explicit vec2d() : x(0), y(0) {}

      vec2d(const T& x_, const T& y_) : x(x_), y(y_) {}

      template <typename S>
            vec2d(const vec2d<S>& o) : x(T(o.x)), y(T(o.y)) {}

      vec2d<T> operator+(const vec2d<T>& p) const {
         return vec2d<T>(x + p.x, y + p.y);
      }

      vec2d<T> operator-(const vec2d<T>& p) const {
         return vec2d<T>(x - p.x, y - p.y);
      }

      template <typename Scalar>
            friend vec2d<T> operator/(const vec2d<T>& p, const Scalar& s) {
         return vec2d<T>(p.x / T(s), p.y / T(s));
      }

      bool operator==(const vec2d& o) const { return x==o.x && y==o.y; }
      bool operator!=(const vec2d& o) const { return !(*this==o); }

      friend std::ostream& operator<<(std::ostream& os, const vec2d<T>& p) {
         return (os << p.x << " " << p.y);
      }

   };

   typedef vec2d<double> point2d;

   template <typename T>
   struct vec3d
   {
#ifdef CVLAB_USE_ARRAY_VEC
      //FIXME: non sembra funzionare
      T xyz[3];
      T &x, &y, &z;

      explicit vec3d() : x(xyz[0]), y(xyz[1]), z(xyz[2]) { x = 0; y = 0; z = 0; }
      vec3d(T x_, T y_, T z_) : x(xyz[0]), y(xyz[1]), z(xyz[2]) { x = x_; y = y_; z = z_; }

      template <typename S>
      vec3d(const vec3d<S>& s) : x(xyz[0]), y(xyz[1]), z(xyz[2]) { x = T(s.x); y = T(s.y); z = T(s.z); }

      const T* to_ptr() const { return xyz; }
      T* to_ptr() { return xyz; }

      vec3d<T>& operator=(const vec3d<T>& p) // needed because of the T&
      {
         x = p.x; y = p.y; z = p.z;
         return *this;
      }
#else
      T x,y,z;

      explicit vec3d() : x(0), y(0), z(0) {}
      vec3d(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

      template <typename S>
      vec3d(const vec3d<S>& s) : x(T(s.x)), y(T(s.y)), z(T(s.z)) {}

      template <typename S>
      vec3d(const S* s) : x(T(s[0])), y(T(s[1])), z(T(s[2])) {}
#endif

      vec3d<T> operator+(const vec3d<T>& p) const
      {
         return vec3d<T>(x + p.x, y + p.y, z + p.z);
      }

      vec3d<T> operator-(const vec3d<T>& p) const
      {
         return vec3d<T>(x - p.x, y - p.y, z - p.z);
      }

      vec3d<T> operator-() const
      {
         return vec3d<T>(-x, -y, -z);
      }

      template <typename S>
      vec3d<T>& operator+=(const vec3d<S>& p) {
         x += T(p.x);
         y += T(p.y);
         z += T(p.z);
         return *this;
      }

      // rules for partial ordering of function templates say that the new overload
      // is a better match, when it matches

      vec3d<T>& operator+=(const vec3d<T>& p)
      {
         x += p.x;
         y += p.y;
         z += p.z;
         return *this;
      }

      template <typename S>
      vec3d<T>& operator-=(const vec3d<S>& p)
      {
         x -= T(p.x);
         y -= T(p.y);
         z -= T(p.z);
         return *this;
      }

      vec3d<T>& operator-=(const vec3d<T>& p) {
         x -= p.x;
         y -= p.y;
         z -= p.z;
         return *this;
      }

      template <typename Scalar>
      vec3d<T>& operator/=(const Scalar& s) {
         const T i = T(1) / T(s);
         x *= i;
         y *= i;
         z *= i;
         return *this;
      }

      template <typename Scalar>
      vec3d<T>& operator*=(const Scalar& s) {
         x *= T(s);
         y *= T(s);
         z *= T(s);
         return *this;
      }

      bool operator==(const vec3d& o) const { return (x == o.x) && (y == o.y) && (z == o.z); }

      template <typename S>
      bool operator==(const vec3d<S>& o) const {
         return (x == T(o.x) && y == T(o.y) && z == T(o.z));
      }

      bool operator!=(const vec3d& o) const { return !(*this==o); }

      template <typename S>
      bool operator!=(const vec3d<S>& o) const {
         return !(*this == o);
      }

      template <typename Scalar>
      friend vec3d<T> operator*(const vec3d<T>& p, const Scalar& s) {
         return vec3d<T>(s * p.x, s * p.y, s * p.z);
      }

      template <typename Scalar>
      friend vec3d<T> operator*(const Scalar& s, const vec3d<T>& p) {
         return p*s;
      }

      template <typename Scalar>
      friend vec3d<T> operator/(const vec3d<T>& p, const Scalar& s) {
         return vec3d<T>(p.x / T(s), p.y / T(s), p.z / T(s));
      }

      friend std::ostream& operator<<(std::ostream& os, const vec3d<T>& p) {
         return (os << p.x << " " << p.y << " " << p.z);
      }

      friend std::istream& operator>>(std::istream& is, vec3d<T>& p) {
         return (is >> p.x >> p.y >> p.z);
      }

   };

   typedef vec3d<double> normal3d;
   typedef vec3d<double> point3d;


   // FIXME: is this really needed? where am i using this?
   template <typename T>
   class edge3d
   {
   private:
      mutable T len; // cache
   public:
      vec3d<T> src;
      vec3d<T> dest;

      explicit edge3d() : len(-1.) {}
      edge3d(const vec3d<T>& s, const vec3d<T>& d) : len(-1.), src(s), dest(d) {}

      T length(bool recompute=false) const {
         if (recompute || (len == -1.))
            len = dist(src,dest);
         return len;
      }
   };


   // FIXME: refactor?
   class oriented_point3d : public point3d {
   public:
      normal3d n;

      explicit oriented_point3d() : point3d() {}
      oriented_point3d(const point3d& p) : point3d(p) {}
      oriented_point3d(double xx, double yy, double zz) : point3d(xx,yy,zz) {}
      oriented_point3d(const oriented_point3d& p) : point3d(p), n(p.n) {}
      oriented_point3d(const point3d& p, const normal3d& nn) : point3d(p), n(nn) {}
   };


   class invalid_vector : public std::logic_error
   {
   public:
      explicit invalid_vector() : std::logic_error("Exception invalid_vector caught.") {}
      invalid_vector(const std::string& msg) : std::logic_error("Exception invalid_vector caught: "+msg) {}
   };


   class collinear_points : public std::runtime_error
   {
   public:
      explicit collinear_points() : std::runtime_error("Exception collinear_points caught.") {}
      collinear_points(const std::string& msg) : std::runtime_error("Exception collinear_points caught: "+msg) {}
   };


   // ==============================================================
   //                         Rotations
   // ==============================================================

   // NOTE: this is a std::vector derivation, thus a matrix<bool> will
   // take 1 bit per element.

   template<class T>
   class matrix : private std::vector<T> // row-major order
   {
   private:
      typedef std::vector<T> super;
      int width_;
      int height_;

   public:

      const int& width;
      const int& height;

      typedef typename super::iterator iterator;
      typedef typename super::const_iterator const_iterator;

      explicit matrix() : super(), width_(0), height_(0), width(width_), height(height_) {}

      matrix(int w, int h) : super(w*h), width_(w), height_(h), width(width_), height(height_) {}

      matrix(int w, int h, const T& v) : super(w*h, v), width_(w), height_(h), width(width_), height(height_) {}

      matrix(const matrix<T>& m) : super(), width(width_), height(height_)
      {
         resize(m.width_,m.height_);
         super::assign( m.begin(), m.end() );
      }

      typename super::reference operator() (size_t r, size_t c) { return super::operator[](r*width_+c); }
      typename super::const_reference operator() (size_t r, size_t c) const { return super::operator[](r*width_+c); }

      typename super::reference at(const size_t r, const size_t c) { return super::at(r*width_+c); }
      typename super::const_reference at(size_t r, size_t c) const { return super::at(r*width_+c); }

      const T* to_ptr() const {
         return &(super::operator[](0)); // ok since std::vector is guaranteed to be contiguous
      }

      T* to_ptr() {
         return &(super::operator[](0));
      }

      void resize(int w, int h)
      {
         super::resize(w*h);
         width_ = w;
         height_ = h;
      }

      size_t size() const { return super::size(); }

      // invalidates all previously obtained iterators, references and pointers
      iterator add_row(T val=0)
      {
         super::insert( std::vector<T>::end(), width_, val );
         return (super::end() - width_); //TESTME
      }

      void clear() { super::clear(); }

      iterator begin() { return super::begin(); }
      const_iterator begin() const { return super::begin(); }
      iterator end() { return super::end(); }
      const_iterator end() const { return super::end(); }

      iterator row(int k) { return (super::begin() + k * width); }
      const_iterator row(int k) const { return (super::begin() + k * width); }

      bool operator==(const matrix<T>& m) const
      {
         if ((width_ != m.width_) || (height != m.height_ ))
            return false;

         const_iterator it1(begin()), it2(m.begin()), it1_end(end());

         for(; it1 != it1_end; ++it1, ++it2)
            if (*it1 != *it2) return false;

         return true;
      }

      bool operator!=(const matrix<T>& m) const { return !(*this == m); }

      matrix& operator=(const matrix<T>& m)
      {
         if (&m == this)
            return *this;

         if (width != m.width || height != m.height)
            throw invalid_vector("Cannot assign matrices with different sizes.");

         super::assign( m.begin(), m.end() );
         return *this;
      }

      template <typename S>
      matrix<T>& operator*=(const S& s)
      {
         for (size_t i=0; i<size(); ++i)
            super::operator[](i) *= T(s);
         return *this;
      }

      template <typename S>
      matrix<T>& operator/=(const S& s)
      {
         for (size_t i=0; i<size(); ++i)
            super::operator[](i) /= T(s);
         return *this;
      }

      friend std::ostream& operator<<(std::ostream& s, const matrix<T>& m)
      {
         for (int y=0; y<m.height_; ++y) {
            for (int x=0; x<m.width_; ++x) {
               s << m[y*m.width_+x] << " ";
            }
            s << std::endl;
         }
         return s;
      }
   };


   template<typename T>
   struct matrix3x3
   {
      T r00, r01, r02,
      r10, r11, r12,
      r20, r21, r22;

      int width, height;

      explicit matrix3x3()
         : r00(0), r01(0), r02(0)
         , r10(0), r11(0), r12(0)
         , r20(0), r21(0), r22(0)
         , width(3), height(3)
      {}

      template <typename S>
      explicit matrix3x3(const S* v)
         : r00(v[0]), r01(v[1]), r02(v[2])
         , r10(v[3]), r11(v[4]), r12(v[5])
         , r20(v[6]), r21(v[7]), r22(v[8])
         , width(3), height(3)
      {}

      void set_column(size_t c, const vec3d<T>& v)
      {
         T x = v.x;
         T y = v.y;
         T z = v.z;

         if (c==0) {
            r00 = x;
            r10 = y;
            r20 = z;
         }
         else if (c==1) {
            r01 = x;
            r11 = y;
            r21 = z;
         }
         else if (c==2) {
            r02 = x;
            r12 = y;
            r22 = z;
         }
         else
            throw std::logic_error("Cannot set column for 3x3 matrix.");
      }

      T& operator() (size_t row, size_t col)
      {
         switch (row) {
         case 0:
            if (col==0) return r00;
            if (col==1) return r01;
            if (col==2) return r02;
            else throw std::out_of_range("Cannot access element in 3x3 matrix");
         case 1:
            if (col==0) return r10;
            if (col==1) return r11;
            if (col==2) return r12;
            else throw std::out_of_range("Cannot access element in 3x3 matrix");
         case 2:
            if (col==0) return r20;
            if (col==1) return r21;
            if (col==2) return r22;
            else throw std::out_of_range("Cannot access element in 3x3 matrix");
         default:
            throw std::out_of_range("Cannot access element in 3x3 matrix");
         }
      }

      const T& operator() (size_t row, size_t col) const
      {
         switch (row) {
         case 0:
            if (col==0) return r00;
            if (col==1) return r01;
            if (col==2) return r02;
            else throw std::out_of_range("Cannot access element in 3x3 matrix");
         case 1:
            if (col==0) return r10;
            if (col==1) return r11;
            if (col==2) return r12;
            else throw std::out_of_range("Cannot access element in 3x3 matrix");
         case 2:
            if (col==0) return r20;
            if (col==1) return r21;
            if (col==2) return r22;
            else throw std::out_of_range("Cannot access element in 3x3 matrix");
         default:
            throw std::out_of_range("Cannot access element in 3x3 matrix");
         }
      }

      friend std::ostream& operator<<(std::ostream& s, const matrix3x3<T>& m) {
         s << m.r00 << " " << m.r01 << " " << m.r02 << std::endl;
         s << m.r10 << " " << m.r11 << " " << m.r12 << std::endl;
         s << m.r20 << " " << m.r21 << " " << m.r22 << std::endl;
         return s;
      }
   };


   template <typename T>
   inline matrix3x3<T> identity3x3()
   {
      matrix3x3<T> m;
      set_identity(m);
      return m;
   }


   template <typename T>
   inline matrix3x3<T> rotation_x(double angle)
   {
      matrix3x3<T> rot;
      rot.r00 = 1;
      rot.r01 = rot.r02 = rot.r10 = rot.r20 = 0;
      rot.r11 = rot.r22 = std::cos(angle);
      rot.r12 = -std::sin(angle);
      rot.r21 = -rot.r12;
      return rot;
   }

   template <typename T>
   inline matrix3x3<T> rotation_y(double angle)
   {
      matrix3x3<T> rot;
      rot.r01 = rot.r10 = rot.r12 = rot.r21 = 0;
      rot.r11 = 1;
      rot.r00 = rot.r22 = std::cos(angle);
      rot.r02 = std::sin(angle);
      rot.r20 = -rot.r02;
      return rot;
   }

   template <typename T>
   inline matrix3x3<T> rotation_z(double angle)
   {
      matrix3x3<T> rot;
      rot.r02 = rot.r12 = rot.r20 = rot.r21 = 0;
      rot.r22 = 1;
      rot.r00 = rot.r11 = std::cos(angle);
      rot.r01 = -std::sin(angle);
      rot.r10 = -rot.r01;
      return rot;
   }


   template <typename T>
   inline void mult_matrix_inplace(const matrix3x3<T>& m1, const matrix3x3<T>& m2, matrix3x3<T>& r)
   {
      const double r00 = m1.r00*m2.r00 + m1.r01*m2.r10 + m1.r02*m2.r20;
      const double r01 = m1.r00*m2.r01 + m1.r01*m2.r11 + m1.r02*m2.r21;
      const double r02 = m1.r00*m2.r02 + m1.r01*m2.r12 + m1.r02*m2.r22;

      const double r10 = m1.r10*m2.r00 + m1.r11*m2.r10 + m1.r12*m2.r20;
      const double r11 = m1.r10*m2.r01 + m1.r11*m2.r11 + m1.r12*m2.r21;
      const double r12 = m1.r10*m2.r02 + m1.r11*m2.r12 + m1.r12*m2.r22;

      const double r20 = m1.r20*m2.r00 + m1.r21*m2.r10 + m1.r22*m2.r20;
      const double r21 = m1.r20*m2.r01 + m1.r21*m2.r11 + m1.r22*m2.r21;
      const double r22 = m1.r20*m2.r02 + m1.r21*m2.r12 + m1.r22*m2.r22;

      r.r00 = r00; r.r01 = r01; r.r02 = r02;
      r.r10 = r10; r.r11 = r11; r.r12 = r12;
      r.r20 = r20; r.r21 = r21; r.r22 = r22;
   }

   template <typename T>
   inline void mult_matrix(const matrix3x3<T>& m1, const matrix3x3<T>& m2, matrix3x3<T>& r)
   {
      if (&r == &m1 || &r==&m2) throw std::logic_error("cvlab::mult_matrix() argument alias");
      return mult_matrix_inplace<T>(m1,m2,r);
   }

   // NO in-place!
   template <typename Rot1, typename Rot2, typename Rot3>
   void mult_matrix(const Rot1& m1, const Rot2& m2, Rot3& r)
   {
      if ((char*)&r == (char*)&m1 || (char*)&r == (char*)&m2)
         throw std::logic_error("cvlab::mult_matrix() argument alias");

      if (m1.width != m2.height || r.height != m1.height || r.width != m2.width)
         throw std::logic_error("Incompatible size matrices");

      double sum;

      for (int is=0; is<m1.height; ++is)
      {
         for (int jd=0; jd<m2.width; ++jd)
         {
            sum = 0.;
            for (int js=0; js<m1.width; ++js)
            {
               sum += m1(is,js) * m2(js,jd);
            }
            r(is,jd) = sum;
         }
      }
   }

   // ==============================================================
   //                         Quaternions
   // ==============================================================

   template <typename T>
   struct quaternion
   {
      T w, i, j, k;

      explicit quaternion(T v=0) : w(v), i(0), j(0), k(0) {}
      //explicit quaternion(const T* p) : w(p[0]), i(p[1]), j(p[2]), k(p[3]) {}
      quaternion(T ww, T ii, T jj, T kk) : w(ww), i(ii), j(jj), k(kk) {}

      static quaternion<T> convert(const vec3d<T>& p) { return quaternion<T>(0, p.x, p.y, p.z); }

      static quaternion<T> convert(const T* p) { return quaternion<T>(p[0], p[1], p[2], p[3]); }

      quaternion<T>& operator+= (const quaternion<T>& a) {w+=a.w; i+=a.i; j+=a.j; k+=a.k; return *this;}

      quaternion<T>& operator*= (T a) {w*=a; i*=a; j*=a; k*=a; return *this;}

      void to_vector(T* p) const { p[0]=w; p[1]=i; p[2]=j; p[3]=k; }

      void to_axis(T* axis) const
      {
         const T n2 = i*i + j*j + k*k;
         if (cvlab::almost_zero(n2,1e-8))
         {
            axis[0] = axis[1] = axis[2] = 0;
            return;
         }

         double mag = std::sqrt(n2);
         mag = (mag == cvlab::almost_zero(mag,1e-8) ? 1. : mag);
         axis[0] = - i / mag;
         axis[1] = - j / mag;
         axis[2] = - k / mag;
      }

      friend std::ostream& operator<<(std::ostream& os, const quaternion<T>& q)
      {
         return os << "[ " << q.w << " " << q.i << " " << q.j << " " << q.k << " ]";
      }

      friend std::istream& operator>>(std::istream& is, quaternion<T>& q)
      {
         std::string dump;
         return (is >> dump >> q.w >> q.i >> q.j >> q.k >> dump);
      }
   };

   template <typename T>
   quaternion<T> operator+ (const quaternion<T>& a, const quaternion<T>& b)
   {
      quaternion<T> result(a);
      result += b;
      return result;
   }

   template <typename T>
   T dot(const quaternion<T>& a, const quaternion<T>& b)
   {
      return a.w*b.w + a.i*b.i + a.j*b.j + a.k*b.k;
   }

   template <typename T>
   T norm(const quaternion<T>& a)
   {
      return std::sqrt(dot(a,a));
   }

   template <typename T>
   quaternion<T> operator* (const quaternion<T>& a, const quaternion<T>& b)
   {
      quaternion<T> result;

      result.w = a.w*b.w - a.i*b.i - a.j*b.j - a.k*b.k;
      result.i = a.i*b.w + a.w*b.i - a.k*b.j + a.j*b.k;
      result.j = a.j*b.w + a.k*b.i + a.w*b.j - a.i*b.k;
      result.k = a.k*b.w - a.j*b.i + a.i*b.j + a.w*b.k;

      return result;
   }

   template <typename T>
   quaternion<T> operator~ (const quaternion<T>& a)
   {
      return quaternion<T>(a.w, -a.i, -a.j, -a.k);
   }

   template <typename T>
   inline void conjugate(quaternion<T>& q)
   {
      q.i = -q.i;
      q.j = -q.j;
      q.k = -q.k;
   }

   template <typename T>
   inline void normalize(quaternion<T>& q)
   {
      T mag = q.w*q.w + q.i*q.i + q.j*q.j + q.k*q.k;
      if (!almost_zero(mag-T(1), 1e-10))
      {
         mag = std::sqrt(mag);
         q.w /= mag;
         q.i /= mag;
         q.j /= mag;
         q.k /= mag;
      }
   }

   template <typename T>
   inline void set_identity(quaternion<T>& q)
   {
      q.w = T(1);
      q.i = q.j = q.k = T(0);
   }

   // returns a normalized unit quaternion
   template<typename T>
   quaternion<T> rot_matrix_to_quaternion(const matrix3x3<T>& m)
   {
      const T m00 = m(0,0);
      const T m11 = m(1,1);
      const T m22 = m(2,2);
      const T m01 = m(0,1);
      const T m02 = m(0,2);
      const T m10 = m(1,0);
      const T m12 = m(1,2);
      const T m20 = m(2,0);
      const T m21 = m(2,1);
      const T tr = m00 + m11 + m22;

      quaternion<T> ret;

      if (tr > 0)
      {
         T s = std::sqrt(tr+T(1)) * 2; // S=4*qw
         ret.w = 0.25 * s;
         ret.i = (m21 - m12) / s;
         ret.j = (m02 - m20) / s;
         ret.k = (m10 - m01) / s;
      }
      else if ((m00 > m11)&(m00 > m22))
      {
         const T s = std::sqrt(T(1) + m00 - m11 - m22) * 2; // S=4*qx
         ret.w = (m21 - m12) / s;
         ret.i = 0.25 * s;
         ret.j = (m01 + m10) / s;
         ret.k = (m02 + m20) / s;
      }
      else if (m11 > m22)
      {
         const T s = std::sqrt(T(1) + m11 - m00 - m22) * 2; // S=4*qy
         ret.w = (m02 - m20) / s;
         ret.i = (m01 + m10) / s;
         ret.j = 0.25 * s;
         ret.k = (m12 + m21) / s;
      }
      else
      {
         const T s = std::sqrt(T(1) + m22 - m00 - m11) * 2; // S=4*qz
         ret.w = (m10 - m01) / s;
         ret.i = (m02 + m20) / s;
         ret.j = (m12 + m21) / s;
         ret.k = 0.25 * s;
      }

      return ret;
   }

   // assumes a normalized unit quaternion
   template <typename T>
   matrix3x3<T> quaternion_to_rot_matrix(const quaternion<T>& q)
   {
      matrix3x3<T> m;
      const T w = q.w, i = q.i, j = q.j, k = q.k;

      m(0,0) = 1 - 2*j*j - 2*k*k;
      m(0,1) = 2*i*j - 2*k*w;
      m(0,2) = 2*i*k + 2*j*w;

      m(1,0) = 2*i*j + 2*k*w;
      m(1,1) = 1 - 2*i*i - 2*k*k;
      m(1,2) = 2*j*k - 2*i*w;

      m(2,0) = 2*i*k - 2*j*w;
      m(2,1) = 2*j*k + 2*i*w;
      m(2,2) = 1 - 2*i*i - 2*j*j;

      return m;
   }

   //TESTME
   template <typename T>
   inline void mult_quaternion(const quaternion<T>& a, const quaternion<T>& b, quaternion<T>& r)
   {
      r.w = a.w*b.w - (a.i*b.i + a.j*b.j + a.k*b.k);
      r.i = a.w*b.i + b.w*a.i + a.j*b.k - a.k*b.j;
      r.j = a.w*b.j + b.w*a.j + a.k*b.i - a.i*b.k;
      r.k = a.w*b.k + b.w*a.k + a.i*b.j - a.j*b.i;
   }

   // ==============================================================
   //
   // ==============================================================

   template <typename T>
   void transpose(matrix<T>& m)
   {
      matrix<T> old(m);
      const int w = m.width;
      const int h = m.height;
      m.resize(h,w);

      for (int row=0; row<h; ++row) {
         for (int col=0; col<w; ++col) {
            m(col,row) = old(row,col);
         }
      }
   }

   template <typename T>
   inline void transpose(matrix3x3<T>& m)
   {
      const T m01 = m.r01, m02 = m.r02, m12 = m.r12, m10 = m.r10, m21 = m.r21, m20 = m.r20;
      m.r01 = m10; m.r02 = m20;
      m.r10 = m01; m.r20 = m02;
      m.r12 = m21; m.r21 = m12;
   }

   template <typename T>
   inline matrix3x3<T> get_transpose(const matrix3x3<T>& m)
   {
      matrix3x3<T> ret;
      ret.r00 = m.r00; ret.r01 = m.r10; ret.r02 = m.r20;
      ret.r10 = m.r01; ret.r11 = m.r11; ret.r20 = m.r02;
      ret.r12 = m.r21; ret.r21 = m.r12; ret.r22 = m.r22;
      return ret;
   }

   // dest matrix must already be of the right size
   template <typename T>
   void transpose(const matrix<T>& src, matrix<T>& dest)
   {
      if (src.width != dest.height || src.height != dest.width)
         throw cvlab::invalid_vector("cvlab::transpose(): Destination matrix must be of the right size.");

      const int w = src.width;
      const int h = src.height;

      for (int row=0; row<h; ++row) {
         for (int col=0; col<w; ++col) {
            dest(col,row) = src(row,col);
         }
      }
   }

   template <typename T>
   inline void transpose(const matrix3x3<T>& src, matrix3x3<T>& dest)
   {
      dest.r00 = src.r00; dest.r11 = src.r11; dest.r22 = src.r22;
      dest.r01 = src.r10; dest.r02 = src.r20;
      dest.r10 = src.r01; dest.r12 = src.r21;
      dest.r21 = src.r12; dest.r20 = src.r02;
   }

   template <typename T>
   void set_identity(matrix<T>& m, T val=1)
   {
      if (m.width != m.height)
         throw invalid_vector("Cannot set identity on a rectangular matrix.");

      if (m.width == 0)
         return;

      const int n = m.width * m.height;
      const int w = m.width;

      int one = 0;
      for (int k=0; k<n; ++k)
      {
         if (k == one) {
            m(k / w, k % w) = val;
            one += w+1;
         }
         else
            m(k / w, k % w) = 0;
      }
   }

   template <typename T>
   void set_identity(matrix3x3<T>& m, T val=1)
   {
      m.r00 = val; m.r01 = 0; m.r02 = 0;
      m.r10 = 0; m.r11 = val; m.r12 = 0;
      m.r20 = 0; m.r21 = 0; m.r22 = val;
   }

   template <typename T>
   inline void set_diag(matrix<T>& m, T val)
   {
      const int w = m.width;
      if (w != m.height)
         throw invalid_vector("Cannot set trace on a rectangular matrix.");

      for (int one=0; one<w; ++one)
         m(one,one) = val;
   }

   template <typename T>
   inline void set_diag(matrix3x3<T>& m, T val)
   {
      m.r00 = val;
      m.r11 = val;
      m.r22 = val;
   }

   template <typename T>
   inline void set_zero(matrix<T>& m)
   {
      std::fill(m.begin(), m.end(), 0);
   }

   template <typename T>
   inline void set_zero(matrix3x3<T>& m)
   {
      m.r00 = 0; m.r01 = 0; m.r02 = 0;
      m.r10 = 0; m.r11 = 0; m.r12 = 0;
      m.r20 = 0; m.r21 = 0; m.r22 = 0;
   }


   // ==============================================================
   //                         Rigid Motion
   // ==============================================================

   template <typename T>
   void normalize_rot_matrix(matrix3x3<T>& m)
   {
      quaternion<T> q( rot_matrix_to_quaternion(m) );
      normalize(q);
      m = quaternion_to_rot_matrix(q);
   }

   template <typename T>
   class rigid_transform
   {
      matrix3x3<T> rot;
      point3d trans;

      public:

         explicit rigid_transform() { set_identity(rot); }

         rigid_transform(const matrix3x3<T>& rot_, const point3d& trans_) : rot(rot_), trans(trans_) {}

         rigid_transform(const matrix<T> rot_, const point3d& trans_) : trans(trans_)
         {
            if (rot_.height != 3 || rot_.width != 3)
               throw std::logic_error("The rotation matrix must be 3x3.");

            rot.r00 = rot_(0,0); rot.r01 = rot_(0,1); rot.r02 = rot_(0,2);
            rot.r10 = rot_(1,0); rot.r11 = rot_(1,1); rot.r12 = rot_(1,2);
            rot.r20 = rot_(2,0); rot.r21 = rot_(2,1); rot.r22 = rot_(2,2);
         }

         void invert()
         {
            transpose(rot);
            rotate(trans, rot);
            trans *= T(-1);
         }

         rigid_transform get_invert() const
         {
            matrix3x3<T> ri = get_transpose(rot);
            point3d ti = get_rotate(trans, ri);
            ti *= T(-1);
            return rigid_transform(ri, ti);
         }

         const matrix3x3<T>& R() const { return rot; }
         const point3d& t() const { return trans; }
   };

   template <typename T>
   void rotate(vec3d<T>& p, const matrix<T>& rot)
   {
      if (rot.height != 3 || rot.width != 3)
         throw std::logic_error("Rotation matrix must be 3x3");

      T oldx = p.x, oldy = p.y, oldz = p.z;
      p.x = oldx*rot(0,0) + oldy*rot(0,1) + oldz*rot(0,2);
      p.y = oldx*rot(1,0) + oldy*rot(1,1) + oldz*rot(1,2);
      p.z = oldx*rot(2,0) + oldy*rot(2,1) + oldz*rot(2,2);
   }

   template <typename T>
   void rotate(vec3d<T>& p, const matrix3x3<T>& rot)
   {
      T oldx = p.x, oldy = p.y, oldz = p.z;
      p.x = oldx*rot.r00 + oldy*rot.r01 + oldz*rot.r02;
      p.y = oldx*rot.r10 + oldy*rot.r11 + oldz*rot.r12;
      p.z = oldx*rot.r20 + oldy*rot.r21 + oldz*rot.r22;
   }

   template <typename T, typename S>
   void rotate(vec3d<T>& p, const matrix<S>& rot)
   {
      if (rot.height != 3 || rot.width != 3)
         throw std::logic_error("Rotation matrix must be 3x3");

      T oldx = p.x, oldy = p.y, oldz = p.z;
      p.x = T( oldx*rot(0,0) + oldy*rot(0,1) + oldz*rot(0,2) );
      p.y = T( oldx*rot(1,0) + oldy*rot(1,1) + oldz*rot(1,2) );
      p.z = T( oldx*rot(2,0) + oldy*rot(2,1) + oldz*rot(2,2) );
   }

   template <typename T, typename S>
   inline void rotate(vec3d<T>& p, const matrix3x3<S>& rot)
   {
      T oldx = p.x, oldy = p.y, oldz = p.z;
      p.x = T( oldx*rot.r00 + oldy*rot.r01 + oldz*rot.r02 );
      p.y = T( oldx*rot.r10 + oldy*rot.r11 + oldz*rot.r12 );
      p.z = T( oldx*rot.r20 + oldy*rot.r21 + oldz*rot.r22 );
   }

   //TODO
   template <typename T>
   inline void rotate(vec3d<T>& p, const quaternion<T>& rot)
   {
      rotate(p, quaternion_to_rot_matrix(rot));
   }


   //TESTME
   template <typename T>
   inline vec3d<T> get_rotate(const vec3d<T>& v, const quaternion<T>& q)
   {
      return get_rotate(v, quaternion_to_rot_matrix(q));
      /*
      const T
            a = q.w, b = q.i, c = q.j, d = q.k,
            t2 = a*b, t3 = a*c, t4 = a*d, t5 = -b*b, t6 = b*c,
            t7 = b*d, t8 = -c*c, t9 = c*d, t10 = -d*d,
            v1 = v.x, v2 = v.y, v3 = v.z;
      return vec3d<T>(
            2*( (t8 + t10)*v1 + (t6 -  t4)*v2 + (t3 + t7)*v3 ) + v1,
            2*( (t4 +  t6)*v1 + (t5 + t10)*v2 + (t9 - t2)*v3 ) + v2,
            2*( (t7 -  t3)*v1 + (t2 +  t9)*v2 + (t5 + t8)*v3 ) + v3);
            */
   }

   template <typename T>
   inline vec3d<T> get_rotate(const vec3d<T>& p, const matrix3x3<T>& rot)
   {
      return vec3d<T>(p.x*rot.r00 + p.y*rot.r01 + p.z*rot.r02,
                      p.x*rot.r10 + p.y*rot.r11 + p.z*rot.r12,
                      p.x*rot.r20 + p.y*rot.r21 + p.z*rot.r22);
   }


   template <typename T, typename RotationType>
   inline void rotate_translate(vec3d<T>& v, const RotationType& rot, const point3d& trans)
   {
      rotate(v, rot);
      v += trans;
   }

   template <typename T>
   inline vec3d<T> get_rotate_translate(const vec3d<T>& p, const matrix3x3<T>& rot, const point3d& t)
   {
      return vec3d<T>(p.x*rot.r00 + p.y*rot.r01 + p.z*rot.r02 + t.x,
                      p.x*rot.r10 + p.y*rot.r11 + p.z*rot.r12 + t.y,
                      p.x*rot.r20 + p.y*rot.r21 + p.z*rot.r22 + t.z);
   }

   template <typename T>
   inline vec3d<T> get_rotate_translate(const vec3d<T>& p, const matrix<T>& rot, const point3d& t)
   {
      if (rot.height != 3 || rot.width != 3)
         throw std::logic_error("Rotation matrix must be 3x3");

      return vec3d<T>(p.x*rot(0,0) + p.y*rot(0,1) + p.z*rot(0,2) + t.x,
                      p.x*rot(1,0) + p.y*rot(1,1) + p.z*rot(1,2) + t.y,
                      p.x*rot(2,0) + p.y*rot(2,1) + p.z*rot(2,2) + t.z);
   }

   template <typename T>
   inline vec3d<T> get_rotate_translate(const vec3d<T>& p, const T* rot, const T* t)
   {
      return get_rotate_translate(p, matrix3x3<T>(rot), vec3d<T>(t[0],t[1],t[2]));
   }

   template <typename T>
   inline vec3d<T> get_rotate_translate(const vec3d<T>& v, const quaternion<T>& rot, const point3d& t)
   {
      return (get_rotate(v, rot) + t);
   }


   /**
    * Inverts a rigid motion.
    */
   template <typename R, typename T>
   inline void invert(R& r, T& t)
   {
      transpose(r);
      rotate(t, r);
      t.x = -t.x; t.y = -t.y; t.z = -t.z;
   }


   /**
    * Computes the rigid motion bringing points expressed in j-coordinates
    * towards the world i, i.e.: Pi = Rij * Pj + Tij
    */
   template <typename T>
   void relative_motion(

         const matrix3x3<T>& Ri, const point3d& Ti,
         const matrix3x3<T>& Rj, const point3d& Tj,
         matrix3x3<T>& Rij, point3d& Tij

   ){
      matrix3x3<T> Ri_inv = Ri;
      point3d Ti_inv = Ti;
      invert(Ri_inv, Ti_inv);

      mult_matrix(Ri_inv, Rj, Rij);
      Tij = get_rotate_translate(Tj, Ri_inv, Ti_inv);
   }


   // ==============================================================
   //                      Vector Operations
   // ==============================================================

   template <typename T>
   inline double normalize(vec3d<T>& p)
   {
      const double n = magnitude(p);
      if (n==0.) {
         p.x = 0;
         p.y = 0;
         p.z = 0;
      }
      else {
         p.x /= n;
         p.y /= n;
         p.z /= n;
      }
      return n;
   }

   template <typename T>
   inline vec3d<T> get_normalize(const vec3d<T>& p)
   {
      vec3d<T> q(p);
      normalize(q);
      return q;
   }

   template <typename T>
   inline double dist(const T& p1, const T& p2)
   {
      const double sqdist = squared_dist(p1,p2);
      return (sqdist == 0. ? 0. : std::sqrt(sqdist));
   }

   // ||p1 - p2||^2
   template <typename T>
         inline double squared_dist(const vec3d<T>& p1, const vec3d<T>& p2)
   {
      T x = p1.x - p2.x;
      T y = p1.y - p2.y;
      T z = p1.z - p2.z;
      return ((x*x) + (y*y) + (z*z));
   }

   template <typename T>
   inline double squared_dist(const vec2d<T>& p1, const vec2d<T>& p2)
   {
      T x = p1.x - p2.x;
      T y = p1.y - p2.y;
      return ((x*x) + (y*y));
   }

   template <typename T>
         inline T dot_product(const vec2d<T>& v1, const vec2d<T>& v2) {
      return (v1.x*v2.x) + (v1.y*v2.y);
   }

   template <typename T>
         inline T dot_product(const vec3d<T>& v1, const vec3d<T>& v2) {
      return (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z);
   }

   template <typename T, typename S>
         inline double dot_product(const vec3d<T>& v1, const vec3d<S>& v2) {
      return (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z);
   }

   template <typename T>
   inline T dot_product(const quaternion<T>& p, const quaternion<T>& q) {
      return (p.w*q.w + p.i*q.i + p.j*q.j + p.k*q.k);
   }

   template <typename T>
   inline double norm2(const T& v)
   {
      return dot_product(v,v);
   }

   template <typename T>
   inline double magnitude(const T& p)
   {
      return std::sqrt(dot_product(p,p));
   }

   template <typename T>
         inline vec3d<T> cross_product(const vec3d<T>& v1, const vec3d<T>& v2)
   {
      return vec3d<T>(
            (v1.y*v2.z) - (v1.z*v2.y),
            (v1.z*v2.x) - (v1.x*v2.z),
            (v1.x*v2.y) - (v1.y*v2.x)
            );
   }


   // ==============================================================
   //                        Misc Geometry
   // ==============================================================

   /**
 * Calculates the absolute distance between two angles.
 * The angles must be given in radians and within [-2pi,2pi].
 */
   template <typename T>
         inline T calc_angle_distance(T alpha, T beta)
   {
      T a = (alpha >= 0. ? alpha : alpha + cvlab::pi2);
      T b = (beta >= 0. ? beta : beta + cvlab::pi2);

      T M, m;

      if (a >= b) {
         M = a;
         m = b;
      }
      else {
         M = b;
         m = a;
      }

      const T Mm = M - m;
      return ( Mm <= cvlab::pi ? Mm : cvlab::pi2 - Mm );
   }


   /**
 * Linear correlation coefficient.
 * Measures the normalized error using the distance between the data and the best
 * least squares fit line to the data.
 *
 * NOTE: The function needs the two containers to be of the same size to give meaningful results.
 *
 * @param first_p
 * @param first_q
 * @param last
 * @return A value between -1 (anti-correlated) and 1 (completely correlated).
 */
   template <typename Iterator>
   double linear_corr_coeff(Iterator first_p, Iterator first_q, Iterator last)
   {
      typedef typename Iterator::value_type vtype;

      vtype sum_pq = 0., sum_pp = 0., sum_qq = 0.;
      uint32_t n = 0;

      // Pearson algorithm ( http://en.wikipedia.org/wiki/Pearson_correlation )
      n = (uint32_t)( last - first_p );
      double mean_x = double(*first_p);
      double mean_y = double(*first_q);
      ++first_p; ++first_q;
      int i = 2;
      for (; first_p != last; ++first_p, ++first_q, ++i) {
         const double sweep = (i - 1.) / i;
         const double delta_x = double(*first_p) - mean_x;
         const double delta_y = double(*first_q) - mean_y;
         sum_pp += delta_x * delta_x * sweep;
         sum_qq += delta_y * delta_y * sweep;
         sum_pq += delta_x * delta_y * sweep;
         mean_x += delta_x / double(i);
         mean_y += delta_y / double(i);
      }
      const double pop_sd_x = std::sqrt(sum_pp / n);
      const double pop_sd_y = std::sqrt(sum_qq / n);
      const double cov = sum_pq / double(n);

      double ret = cov / (pop_sd_x * pop_sd_y);
      if (ret > 1.) ret = 1.;
      if (ret < -1.) ret = -1.;

      return ret;
   }

   template <typename Iterator>
         double linear_corr_coeff(Iterator first_p, Iterator first_q, Iterator last, uint32_t& n)
         //TESTME!!!!!!!!!!
   {
      typedef typename Iterator::value_type vtype;

      vtype sum_pq = 0., sum_pp = 0., sum_qq = 0.;
      n = 0;

      // Pearson algorithm ( http://en.wikipedia.org/wiki/Pearson_correlation )

      double mean_x = double(*first_p);
      double mean_y = double(*first_q);
      ++first_p; ++first_q;
      int i = 2;
      for (; first_p != last; ++first_p, ++first_q, ++i) { //FIXME: should increment "i" always?
         if (*first_p == 0 || *first_q == 0)
            continue;

         const double sweep = (i - 1.) / i;
         const double delta_x = double(*first_p) - mean_x;
         const double delta_y = double(*first_q) - mean_y;
         sum_pp += delta_x * delta_x * sweep;
         sum_qq += delta_y * delta_y * sweep;
         sum_pq += delta_x * delta_y * sweep;
         mean_x += delta_x / double(i);
         mean_y += delta_y / double(i);
         ++n;
      }

      if (n == 0)
         return -1.;

      const double pop_sd_x = std::sqrt(sum_pp / n);
      const double pop_sd_y = std::sqrt(sum_qq / n);
      const double cov = sum_pq / double(n);

      double ret = cov / (pop_sd_x * pop_sd_y);
      if (ret > 1.) ret = 1.;
      if (ret < -1.) ret = -1.;

      return ret;
   }


   // Heron's numerically-stable formula for the area of a triangle
   template <typename T>
   double get_triangle_area(const vec3d<T>& v0, const vec3d<T>& v1, const vec3d<T>& v2)
   {
      const double ta = dist(v0,v1);
      const double tb = dist(v0,v2);
      const double tc = dist(v1,v2);

      double a, b, c; // a >= b >= c
      if (ta >= tb) {
         if (ta >= tc) {
            a = ta;
            if (tb >= tc) { b = tb; c = tc; }
            else { b = tc; c = tb; }
         }
         else { a = tc; b = ta; c = tb; }
      }
      else {
         if (ta >= tc) { a = tb; b = ta; c = tc; }
         else {
            c = ta;
            if (tb >= tc) { a = tb; b = tc; }
            else { a = tc; b = tb; }
         }
      }

      return ( 0.25 * std::sqrt( (a+(b+c)) * (c-(a-b)) * (c+(a-b)) * (a+(b-c)) ) );
   }


   /**
 *     v1     v2
 *      |   /
 *      |  /
 *      |_/  calculates this angle
 *      |/
 *      v0
 */
   template <typename T>
         double get_angle(const vec3d<T>& v1, const vec3d<T>& v0, const vec3d<T>& v2)
   {
      vec3d<T> v21 = v1 - v0;
      vec3d<T> v23 = v2 - v0;

      double v21_l2 = dot_product(v21,v21);
      double v23_l2 = dot_product(v23,v23);

      if( v21_l2 != 0.0 )
         v21 /= std::sqrt(v21_l2);
      if( v23_l2 != 0.0 )
         v23 /= std::sqrt(v23_l2);

      const double bound = 0.999999;
      double cos_theta = std::max( -bound, std::min( bound, dot_product(v21,v23) ) );

      return std::acos(cos_theta);
   }


   template <typename T>
         double cotangent(const vec3d<T>& v1, const vec3d<T>& v2, const vec3d<T>& v3)
   {
      vec3d<T> v21 = v1 - v2;
      double v21_l2 = dot_product(v21,v21);
      if( v21_l2 != 0.0 )
         v21 /= std::sqrt( v21_l2 );
      else {}

      vec3d<T> v23 = v3 - v2;
      double v23_l2 = dot_product(v23,v23);
      if( v23_l2 != 0.0 )
         v23 /= std::sqrt( v23_l2 );
      else {}

      const double bound = 0.999999;
      double cos_theta = std::max(-bound, std::min(bound, dot_product(v21,v23)));

      return 1.0 / std::tan( std::acos(cos_theta) );
   }


   template <typename T>
         void get_barycentric_coords(

               const vec3d<T>& p0, const vec3d<T>& p1, const vec3d<T>& p2,
               const vec3d<T>& px,
               vec3d<T>& coords

               ) {
      vec3d<T> s = cross_product(p1-p0,p2-p0);
      vec3d<T> n = get_normalize(s);

      // Compute twice area of triangle ABC
      const double AreaABC = dot_product(n,s);

      const double AreaPBC = dot_product(n,cross_product(p1-px, p2-px));
      coords.x = T(AreaPBC / AreaABC);

      const double AreaPCA = dot_product(n,cross_product(p2-px, p0-px));
      coords.y = T(AreaPCA / AreaABC);

      coords.z = 1. - coords.x - coords.y;
   }


   template <typename T>
         inline bool collinear(const vec3d<T>& p, const vec3d<T>& q, const vec3d<T>& r)
   {
      double d = fabs(dot_product(get_normalize(q-p), get_normalize(r-q)));
      if (d>1.) d=1.;

      return almost_zero(1.-d , 1e-15);
   }


   template <typename Iterator>
   inline double median(Iterator start, Iterator end)
   {
      const typename Iterator::difference_type n = end - start;
      if (n <= 0) return 0.;

      if (n % 2 == 0)
         return (*(start+(n/2)) + *(start+(n/2-1))) / 2.;
      else
         return *(start + ( (n-1) / 2 ));
   }


   template <typename Iterator>
   inline double mean(Iterator start, Iterator end)
   {
      const typename Iterator::difference_type n = end - start;
      if (n <= 0) return 0.;

      double acc = 0.;

      for(; start != end; ++start)
         acc += double(*start);

      acc /= n;
      return acc;
   }


   //TESTME
   template <typename Iterator>
   inline double variance(Iterator start, Iterator end, double mean)
   {
      double var = 0.;

      for(; start != end; ++start)
         var += (double(*start)-mean) * (double(*start)-mean);

      return var / double(end - start - 1);
   }


   void interpolate_normals(oriented_point3d& p, const oriented_point3d& q0, const oriented_point3d& q1, const oriented_point3d& q2);
   double rand_gauss(double mean, double stdev);
   bool isnan(double var);

   int GCD(int a, int b);
   int LCM(int a, int b);

   bool findLineTriangleIntersection(
         const point3d& S,
         const point3d& V,
         const point3d& P0,
         const point3d& P1,
         const point3d& P2,
         point3d& intersection);

   bool findLineTriangleIntersection2(
         const point3d& S,
         const point3d& V,
         const point3d& P0,
         const point3d& P1,
         const point3d& P2,
         point3d& intersection);

   bool findLineTriangleIntersection3(
         const point3d& S,
         const point3d& V,
         const point3d& P0,
         const point3d& P1,
         const point3d& P2,
         point3d& intersection);

   bool findLineTriangleIntersection4(
         const point3d& S,
         const point3d& V,
         const point3d& P0,
         const point3d& P1,
         const point3d& P2,
         point3d& intersection);

   double findSquareInParallelogram(
         float *A, float *B,
         float *C, float *D);

   point2d IntersectLines(
         const point2d& p1,
         const point2d& p2,
         const point2d& c1,
         const point2d& c2);

   double IntersectLines(
         const point3d& p1,
         const point3d& p2,
         const point3d& c1,
         const point3d& c2,
         double& X, double& Y, double& Z);

   point2d IntersectSegments(
         const point2d& p1,
         const point2d& p2,
         const point2d& p3,
         const point2d& p4);


   /**
    * Provides basic functionalities for pointer vectors.
    */
   namespace ptrmath
   {

      template <typename T>
            inline T dot_product(const T* v1, const T* v2, size_t n)
      {
         T res(0);
         for (size_t i=0; i<n; ++i)
            res += v1[i] * v2[i];
         return res;
      }

      template <typename T>
            inline void cross_product(const T* v1, const T* v2, T* r)
      {
         r[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
         r[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
         r[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
      }

      template <typename T>
            inline T magnitude(const T* v, size_t n)
      {
         return std::sqrt( dot_product(v, v, n) );
      }

      template <typename T>
            inline T norm2(const T* v, size_t n)
      {
         return dot_product(v, v, n);
      }

      // assumes <w i j k>
      template <typename T>
            inline void mult_quaternion(const T* a, const T* b, T* r)
      {
         r[0] = a[0]*b[0] - (a[1]*b[1] + a[2]*b[2] + a[3]*b[3]);
         r[1] = a[0]*b[1] + b[0]*a[1] + a[2]*b[3] - a[3]*b[2];
         r[2] = a[0]*b[2] + b[0]*a[2] + a[3]*b[1] - a[1]*b[3];
         r[3] = a[0]*b[3] + b[0]*a[3] + a[1]*b[2] - a[2]*b[1];
      }

      template <typename T>
            inline void axis_from_quaternion(const T* q, T* axis)
      {
         const T mag = std::sqrt( q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
         axis[0] = - q[1] / mag;
         axis[1] = - q[2] / mag;
         axis[2] = - q[3] / mag;
      }

      template <typename T>
            inline void set_identity(T* mat, size_t n)
      {
         for (size_t i=0; i<n; ++i) {
            for (size_t j=0; j<n; ++j) {
               if (i == j) mat[i*n + j] = T(1);
               else mat[i*n + j] = T(0);
            }
         }
      }

      // does NOT support in-place transposition
      template <typename T>
            void transpose(T *A, T *At, size_t m, size_t n) {
         for (size_t i=0; i<m; ++i)
            for (size_t j=0; j<n; ++j)
               At[j*m + i] = A[i*n + j];
      }

      template <typename T>
            void rotate(T* v, const T* r)
      {
         const T a = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
         const T b = r[3] * v[0] + r[4] * v[1] + r[5] * v[2];
         const T c = r[6] * v[0] + r[7] * v[1] + r[8] * v[2];
         v[0] = a; v[1] = b; v[2] = c;
      }

      template <typename T>
            void rotate_translate(T* v, const T* r, const T* t)
      {
         rotate(v,r);
         v[0] += t[0];
         v[1] += t[1];
         v[2] += t[2];
      }

   }

} // namespace scanner_math


#endif // CVLAB_MATH3D_H
