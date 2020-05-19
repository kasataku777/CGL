#include<iostream>
#include<stdio.h>
#include<math.h>
#include<algorithm>
#include<string>
#include<map>

using namespace std;

#define EPS (1e-10)
#define equals(a,b) (fabs((a)-(b))<EPS)

class Point {
public:
	double x, y;


	Point(double x = 0, double y = 0) :x(x), y(y) {}

	Point operator + (Point p) { return Point(x + p.x, y + p.y); }
	Point operator - (Point p) { return Point(x - p.x, y - p.y); }
	Point operator * (double a) { return Point(a * x, a * y); }
	Point operator / (double a) { return Point(x / a, y / a); }

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
	return equals(dot(a, b), 0.0);
}

bool isParallel(Vector a, Vector b) {
	return equals(cross(a, b), 0.0);
}
struct Segment {
	Point p1, p2;
	Segment() {}
	Segment(Point p1, Point p2) :p1(p1), p2(p2) {}
};
typedef Segment  Line;

double norm(Vector a) {
	return a.x * a.x + a.y * a.y;
}
double abs(Vector a) {
	return sqrt(norm(a));
}

Point project(Segment s, Point p) {
	Vector base = s.p2 - s.p1;
	double r = dot(p - s.p1, base) / norm(base);
	return s.p1 + base * r;
}

Point reflect(Segment s, Point p) {
	return p + (project(s, p) - p) * 2.0;
}
double getDistance(Point a, Point b) {
	return abs(a - b);
}

double getDistanceLP(Line l, Point p) {
	return fabs(cross(l.p2 - l.p1, p - l.p1)/abs(l.p2 - l.p1));
}

double getDistanceSP(Segment s, Point p) {
	if (dot(s.p2 - s.p1, p - s.p1) < 0.0)return abs(p - s.p1);
	if (dot(s.p1 - s.p2, p - s.p2) < 0.0)return abs(p - s.p2);
	return getDistanceLP(s, p);

}


static const int COUNTER_CLOCKWISE = 1;
static const int CLOCKWISE = -1;
static const int ONLINE_BACK = 2;
static const int ONLINE_FRONT = -2;
static const int ON_SEGMENT = 0;

int ccw(Point p0, Point p1, Point p2) {
	Vector a = p1 - p0;
	Vector b = p2 - p0;
	if (cross(a, b) > EPS)return COUNTER_CLOCKWISE;
	if (cross(a, b) < -EPS) return CLOCKWISE;
	if (dot(a, b) < -EPS)return ONLINE_BACK;
	if (a.norm()<b.norm()) return ONLINE_FRONT;

	return ON_SEGMENT;


}
bool intersect(Point p1, Point p2, Point p3, Point p4) {
	return (ccw(p1, p2, p3) * ccw(p1, p2, p4) <= 0 && ccw(p3, p4, p1) * ccw(p3, p4, p2) <= 0);
}
bool intersect(Segment s1, Segment s2) {
	return intersect(s1.p1, s1.p2, s2.p1, s2.p2);
}

double getDistance(Segment s1, Segment s2) {
	if (intersect(s1, s2))return 0.0;
	return min(min(getDistanceSP(s1, s2.p1), getDistanceSP(s1, s2.p2)),
		min(getDistanceSP(s2, s1.p1), getDistanceSP(s2, s1.p2)));
}

Point getCrossPoint(Segment s1, Segment s2) {
	Vector base = s2.p2 - s2.p1;
	double d1 = fabs(cross(base, s1.p1 - s2.p1));
	double d2 = fabs(cross(base, s1.p2 - s2.p1));
	double t = d1 / (d1 + d2);
	return s1.p1 + (s1.p2 - s1.p1) * t;
}

//2B
//int main() {
//	int q,x1,y1,x2,y2,x3,y3,x0,y0;
//	cin >> q;
//	for (int i = 0; i < q; i++) {
//		cin >> x0>>y0>>x1 >> y1 >> x2 >> y2 >> x3 >> y3;
//
//		if (intersect(Point(x0, y0), Point(x1, y1), Point(x2, y2), Point(x3, y3))) {
//			cout << "1" << endl;
//		}
//		else {
//			cout << "0" << endl;
//		}
//	}
//	return 0;
//}

//2D
//int main() {
//	int q;
//	double x1,y1,x2,y2,x3,y3,x0,y0;
//	cin >> q;
//	for (int i = 0; i < q; i++) {
//			cin >> x0>>y0>>x1 >> y1 >> x2 >> y2 >> x3 >> y3;
//			double ans = getDistance(Segment(Point(x0, y0), Point(x1, y1)), Segment(Point(x2, y2), Point(x3, y3)));
//			printf("%.8lf\n", ans);
//	}
//}

//2C

int main() {
	int q;
	cin >> q;
	double x0, y0, x1, y1, x2, y2, x3, y3;
	for (int i = 0; i < q; i++) {
		cin >> x0 >> y0 >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;
		Point ans = getCrossPoint(Segment(Point(x0, y0), Point(x1, y1)), Segment(Point(x2, y2), Point(x3, y3)));
		printf("%.8lf %.8lf\n", ans.x, ans.y);
	}
	return 0;
}