#include<iostream>
#include<math.h>

using namespace std;

#define EPS (1e-10)
#define equals(a,b) (fabs((a)-(b))<EPS)

class Point {
public:
	double x, y;


	Point(double x = 0,double y=0):x(x),y(y){}

	Point operator + (Point p) { return Point(x + p.x, y + p.y); }
	Point operator - (Point p) { return Point(x- p.x, y - p.y); }
	Point operator * (double a) { return Point(a*x, y + a*y); }
	Point operator / (double a) { return Point(x/a, y /a); }

	double abs() { return sqrt(norm()); }
	double norm() { return x * x + y * y; };

	bool operator < (const Point& p)const {
		return x != p.x ? x < p.x : y < p.y;
	}

	bool operator ==(const Point& p)const {
		return fabs(x - p.x) < EPS && fabs(y - p.y) < EPS;
	}

};

typedef Point Vector;

double dot(Vector a, Vector b) {
	return a.x * b.x + a.y * b.y;
}

double cross(Vector a, Vector b) {
	return a.x * b.y - a.y * b.x;
}

bool isOrthogonal(Vector a, Vector b) {
	return equals(dot(a, b) , 0.0);
}

bool isParallel(Vector a, Vector b) {
	return equals(cross(a, b), 0.0);
}

int main(){
	int q;
cin >> q;


for (int i = 0; i < q; i++) {
	int x0, x1, x2, x3, y0, y1, y2, y3;
	cin >> x0 >> y0 >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;

	Vector v0 = Vector(x0, y0);
	Vector v1 = Vector(x1, y1);
	Vector v2 = Vector(x2, y2);
	Vector v3 = Vector(x3, y3);

	Vector s1 = v1 - v0;
	Vector s2 = v3 - v2;

	if (isOrthogonal(s1, s2)) {
		cout << "1" << endl;
	}
	else if (isParallel(s1, s2)) {
		cout << "2" << endl;
	}
	else {
		cout << "0" << endl;
	}

 }

return 0;
}