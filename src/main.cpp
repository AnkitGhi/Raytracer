#include <iostream>
#include <tuple>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <array>    
#include <limits>

using namespace std;
using float4 = std::tuple<float, float, float, float>;
using float3 = std::tuple<float, float, float>;

double roundTo2Decimals(double val) {
    return std::round(val * 100.0) / 100.0;
}

class Tuples{
public:
    float x, y, z, w;

    Tuples(float x_, float y_, float z_, float w_): x(x_),y(y_),z(z_),w(w_) {};

    int point_or_vec() const {
        if (w== 1) cout<<"It's a point. "<<endl;
        else cout<<"It's a vector"<<endl;
        return 0;
    }

    void display() const {
        cout << fixed << setprecision(2);
        cout<<"("<<x<<","<<y<<","<<z<<","<<w<<")"<<endl;
    }

    float magnitude() const {
        float magnitud = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        return magnitud;
    }

    Tuples normalize() const {
        float mag = magnitude();
        return Tuples(x/mag,y/mag,z/mag,w);
    }

    Tuples operator+(const Tuples& second) const{
        return Tuples(x+second.x,y+second.y,z+second.z,w+second.w);
    }

    Tuples operator-(const Tuples& second) const{
        return Tuples(x-second.x,y-second.y,z-second.z,w-second.w);
    }

    Tuples operator-() const{
        return Tuples(-x,-y,-z,-w);
    }

    float operator*(const Tuples& ab) const{
        float dot = x*ab.x + y*ab.y + z*ab.z;
        return dot;
    }

    Tuples operator*(const float num) const{
        return Tuples(x*num,y*num,z*num,w*num);
    }

    Tuples operatorX(const Tuples& ab)const {
        if (w != 0.0f || ab.w != 0.0f){
            throw logic_error("Cross product is only between two vectors");
        }

        return Tuples(
            y * ab.z - z * ab.y,
            z * ab.x - x * ab.z,
            x * ab.y - y * ab.x,
            0.0f
        );
    }
};

Tuples point(float x, float y, float z){
    return Tuples(x, y, z, 1.0);
}

Tuples vector_dec(float x, float y, float z){
    return Tuples(x, y, z, 0.0);
}

int tuple_equality(float3 first, float3 second){
    float epsilon = 1e-6;
    float x = get<0>(first) - get<0>(second);
    float y = get<1>(first) - get<1>(second);
    float z = get<2>(first) - get<2>(second);

    return(x< epsilon && y< epsilon && z< epsilon);
}

float4 add_tuple(float4 first, float4 second){
    float x = get<0>(first) + get<0>(second);
    float y = get<1>(first) + get<1>(second);
    float z = get<2>(first) + get<2>(second);
    float w = get<3>(first) + get<3>(second);
    return make_tuple(x, y, z, w);
}

float4 subtract_tuple(float4 first, float4 second){
    float x = get<0>(first) - get<0>(second);
    float y = get<1>(first) - get<1>(second);
    float z = get<2>(first) - get<2>(second);
    float w = get<3>(first) - get<3>(second);
    return make_tuple(x, y, z, w);
}

float4 negate_tuple(float4 tups){
    float x = 0 - get<0>(tups);
    float y = 0 - get<1>(tups);
    float z = 0 - get<2>(tups);
    float w = 0 - get<3>(tups);
    return make_tuple(x, y, z, w);
}

class projectile{
public:
    Tuples position, velocity;

    projectile(Tuples _position, Tuples _velocity):
        position(_position), velocity(_velocity) {
            if (position.w!=1 && velocity.w!=0){
                throw invalid_argument("Position must be a point and velocity must be a vector");
            }
        };
};

class environment{
public:
    Tuples gravity, wind;

    environment(Tuples& _gravity, Tuples& _wind):
        gravity(_gravity), wind(_wind) {
            if (gravity.w!=0 && wind.w!=0){
                throw invalid_argument("Gravity and wind must be vectors");
            }
        };
};

projectile tick(const projectile& proj, const environment& env){
    Tuples new_position = proj.position + proj.velocity;
    Tuples new_velocity = proj.velocity + env.gravity + env.wind;
    cout<<"The position is ";
    new_position.display();
    return projectile(new_position, new_velocity);
}

class Pixels{
public:
    float red, green, blue;

    Pixels(const float _red, const float _green, const float _blue): 
    red(_red), green(_green), blue(_blue){};

    Pixels operator+(const Pixels& c2)const {
        return Pixels(red+c2.red,green+c2.green, blue+c2.blue);
    }

    Pixels operator-(const Pixels& c2)const {
        return Pixels(red-c2.red,green-c2.green, blue-c2.blue);
    }

    Pixels operator*(const float x)const {
        return Pixels(red*x,green*x, blue*x);
    }

    Pixels operator*(const Pixels& c2)const {
        return Pixels(red*c2.red,green*c2.green, blue*c2.blue);
    }

};

class Canvas{
public:
    int width, height;
    vector<vector <Pixels>> canva_pixels;

    Canvas(int _width, int _height):
    width(_width), height(_height), canva_pixels(height,vector<Pixels>(width,Pixels(0,0,0))){
    };

    void write_pixels(int x, int y, const Pixels& color){
        if (y >= 0 && y < height && x >= 0 && x < width) {
        canva_pixels[y][x] = color;
    }
    }

    string ppm_data(){
        std::stringstream ppm;
        ppm<<"P3\n";
        ppm<<width<<" "<<height<<"\n";
        ppm<<"255\n";

        for(int y = 0; y<height;y++){
            for(int x = 0;x<width;x++){
                int r = clamp(static_cast<int>(round(canva_pixels[y][x].red * 255)),0,255);
                int g = clamp(static_cast<int>(round(canva_pixels[y][x].green * 255)),0,255);
                int b = clamp(static_cast<int>(round(canva_pixels[y][x].blue * 255)),0,255);

                ppm<<r<<" "<<g<<" "<<b<<" ";
            }
            ppm<<"\n";
        }
        return ppm.str();
    }
};

class twomatrix{
public:
    twomatrix(const array<array<float,2>,2>& values): data(values){};

    twomatrix():data({{
        {{1.0f,0.0f}},
        {{0.0f,1.0f}}
    }}){};

    array<float,2>& operator[](int row){
        return data[row];
    }

    const array<float,2>& operator[](int row) const {
        return data[row];
    }

    void print(){
        for (auto &row: data){
            for(auto val: row){
                cout<<val<<" ";
            }
            cout<<endl;
        }
    }

    twomatrix tranpose_matrix()const {
        twomatrix transposed;
        for(int row=0;row<2;row++){
            for(int col=0;col<2;col++){
                transposed[row][col]=data[col][row];
            }
        }
        return transposed;
    }

    twomatrix operator/(float a)const{
        if(a==0){
            throw runtime_error("Cannot divide by zero.");
        }
        twomatrix divided;
        for(int row=0;row<2;row++){
            for(int col=0;col<2;col++){
                divided[row][col]=data[row][col]/a;
            }
        }
        return divided;
    }

    float determinant()const {
        float deter = data[0][0]*(data[1][1])-(data[0][1])*(data[1][0]);
        return deter;
    }

private:
    array<array<float,2>,2> data;
};

class threematrix{
public:
    threematrix(const array<array<float,3>,3>& values): data(values){};
    threematrix(): data({{
        {{1.0f,0.0f,0.0f}},
        {{0.0f,1.0f,0.0f}},
        {{0.0f,0.0f,1.0f}}
    }}){};

    void print(){
        for (auto &row: data){
            for(auto val: row){
                cout<<val<<" ";
            }
            cout<<endl;
        }
    }

    array<float,3>& operator[] (int row){
        return data[row];
    }

    const array<float,3>& operator[] (int row) const{
        return data[row];
    }
    
    threematrix tranpose_matrix()const {
        threematrix transposed;
        for(int row=0;row<3;row++){
            for(int col=0;col<3;col++){
                transposed[row][col]=data[col][row];
            }
        }
        return transposed;
    }

    threematrix operator/(float a)const{
        if(a==0){
            throw runtime_error("Cannot divide by zero.");
        }
        threematrix divided;
        for(int row=0;row<3;row++){
            for(int col=0;col<3;col++){
                divided[row][col]=data[row][col]/a;
            }
        }
        return divided;
    }

    threematrix cofactors()const {
        threematrix cofact;
        for(int x=0;x<3;x++){
            for(int y=0;y<3;y++){
                twomatrix twoma = submatrix(x,y);
                cofact[x][y]=pow(-1,x+y)*twoma.determinant();
            }
        }
        return cofact;
    }

    twomatrix submatrix(int row, int col)const {
        twomatrix newer;
        int a = 0;
        if(row>=3 || col>=3){
            throw out_of_range("Index of the matrices are out of bound. ");
        }
        for(int x = 0;x<3;x++){
            if (x==row){
                ;
            }
            else{
                int b=0;
                for(int y = 0;y<3;y++){
                    if(y==col){
                        ;
                    }
                    else{
                        newer[a][b]=data[x][y];
                        b++;
                    }
                }
                a++;
            }
        }
        return newer;
    }

    float determinant()const {
        twomatrix a = submatrix(0,0);
        twomatrix b = submatrix(0,1);
        twomatrix c = submatrix(0,2);
        float first = data[0][0] * a.determinant();
        float second = data[0][1] * b.determinant();
        float third = data[0][2] * c.determinant();
        float deter = first - second + third;
        return deter;
    }

    threematrix invertible()const {
        float deter = this->determinant();
        if(deter==0){
            throw runtime_error("It is not invertible");
        }

        threematrix cofactor_matrix = this->cofactors();
        threematrix transposed_matrix = cofactor_matrix.tranpose_matrix();
        threematrix inverts = transposed_matrix/deter;
        return inverts;
    }

private:
    array<array<float,3>,3> data;
};

class fourmatrix{
public:
    fourmatrix(const array<array<float,4>,4> &values): data(values){};
    fourmatrix(): data({{
        {{1.0f,0.0f,0.0f,0.0f}},
        {{0.0f,1.0f,0.0f,0.0f}},
        {{0.0f,0.0f,1.0f,0.0f}},
        {{0.0f,0.0f,0.0f,1.0f}}
    }}){};
    
    void print(){
        for (auto &row: data){
            for(auto val: row){
                cout<<val<<" ";
            }
            cout<<endl;
        }
    }

    array<float,4>& operator[](int row){
        return data[row];
    }

    const array<float,4>& operator[](int row)const {
        return data[row];
    }

    int size(){
        return data.size();
    }

    fourmatrix operator*(const fourmatrix& b)const {
        fourmatrix news;
        for(int row=0;row<4;row++){
            for(int col=0;col<4;col++){
                news[row][col]= data[row][0]*b[0][col]+ data[row][1]*b[1][col] + data[row][2]*b[2][col] + data[row][3]*b[3][col];
            }
        }
        return news;
    }

    fourmatrix operator/(float a)const {
        if(a==0){
            throw runtime_error("Cannot divide by zero.");
        }
        fourmatrix divided;
        for(int row=0;row<4;row++){
            for(int col=0;col<4;col++){
                divided[row][col]=data[row][col]/a;
            }
        }
        return divided;
    }

    Tuples operator*(const Tuples& ab)const {
        float a = data[0][0]*ab.x + data[0][1]*ab.y + data[0][2]*ab.z + data[0][3]*ab.w;
        float b = data[1][0]*ab.x + data[1][1]*ab.y + data[1][2]*ab.z + data[1][3]*ab.w;
        float c = data[2][0]*ab.x + data[2][1]*ab.y + data[2][2]*ab.z + data[2][3]*ab.w;
        float d = data[3][0]*ab.x + data[3][1]*ab.y + data[3][2]*ab.z + data[3][3]*ab.w;

        return Tuples(a,b,c,d);
    }

    fourmatrix tranpose_matrix()const{
        fourmatrix transposed;
        for(int row=0;row<4;row++){
            for(int col=0;col<4;col++){
                transposed[row][col]=data[col][row];
            }
        }
        return transposed;
    }

    threematrix submatrix(int row, int col)const {
        threematrix newer;
        int a = 0;
        if(row>=4 || col>=4){
            throw out_of_range("Index of the matrices are out of bound. ");
        }
        for(int x = 0;x<4;x++){
            if (x==row){
                ;
            }
            else{
                int b=0;
                for(int y = 0;y<4;y++){
                    if(y==col){
                        ;
                    }
                    else{
                        newer[a][b]=data[x][y];
                        b++;
                    }
                }
                a++;
            }
        }
        return newer;
    }

    fourmatrix cofactors()const {
        fourmatrix cofact;
        for(int x=0;x<4;x++){
            for(int y=0;y<4;y++){
                threematrix threema = submatrix(x,y);
                cofact[x][y]=pow(-1,x+y)*threema.determinant();
            }
        }
        return cofact;
    }

    float determinant() const {
        threematrix a = submatrix(0,0);
        threematrix b = submatrix(0,1);
        threematrix c = submatrix(0,2);
        threematrix d = submatrix(0,3);
        float first = data[0][0]* a.determinant();
        float second = data[0][1]* b.determinant();
        float third = data[0][2]* c.determinant();
        float fourth = data[0][3]* d.determinant();
        float deter = first - second + third - fourth;
        return deter;
    }

    fourmatrix invertible()const {
        float deter = this->determinant();
        if(deter==0){
            throw runtime_error("It is not invertible");
        }

        fourmatrix cofactor_matrix = this->cofactors();
        fourmatrix transposed_matrix = cofactor_matrix.tranpose_matrix();
        fourmatrix inverts = transposed_matrix/deter;
        return inverts;
    }

private:
    array<array<float,4>,4> data;
};

int matrix_equality(fourmatrix& ab, fourmatrix& bc){
    for (int x= 0;x<4;x++){
        for (int y = 0;y<4;y++){
            if (ab[x][y] != bc[x][y]){
                return 0;
            }
        }
    }
    return 1;
}

fourmatrix translation(float x, float y, float z){
    fourmatrix translate;
    translate[0][0]=1;
    translate[1][1]=1;
    translate[2][2]=1;
    translate[3][3]=1;
    translate[0][3]=x;
    translate[1][3]=y;
    translate[2][3]=z;
    return translate;
}

fourmatrix scaling(float x, float y, float z){
    fourmatrix translate;
    translate[0][0]=x;
    translate[1][1]=y;
    translate[2][2]=z;
    translate[3][3]=1;
    return translate;
}

fourmatrix rotation_x(double angle){
    fourmatrix rota_x;
    rota_x[0][0]=1;
    rota_x[3][3]=1;
    rota_x[1][1]=cos(angle);
    rota_x[2][2]=cos(angle);
    rota_x[2][1]=sin(angle);
    rota_x[1][2]=-sin(angle);
    return rota_x;
}

fourmatrix rotation_y(double angle){
    fourmatrix rota_y;
    rota_y[0][0]=cos(angle);
    rota_y[1][1]=1;
    rota_y[2][2]=cos(angle);
    rota_y[3][3]=1;
    rota_y[2][0]=-sin(angle);
    rota_y[0][2]=sin(angle);
    return rota_y;
}

fourmatrix rotation_z(double angle){
    fourmatrix rota_z;
    rota_z[0][0]=cos(angle);
    rota_z[1][1]=cos(angle);
    rota_z[2][2]=1;
    rota_z[3][3]=1;
    rota_z[0][1]=-sin(angle);
    rota_z[1][0]=sin(angle);
    return rota_z;
}

fourmatrix shearing(int xy, int xz, int yx, int yz, int zx, int zy){
    fourmatrix shear;
    shear[0][0]=1;
    shear[1][1]=1;
    shear[2][2]=1;
    shear[3][3]=1;
    shear[0][1]=xy;
    shear[0][2]=xz;
    shear[1][0]=yx;
    shear[1][2]=yz;
    shear[2][0]=zx;
    shear[2][1]=zy;
    return shear;
}

class Sphere{
public:
    int radius;
    fourmatrix tranform;
    Sphere(int _id): radius(_id),tranform(){}; 

    void set_transform(const fourmatrix& A){
        tranform = A;
    }
};

class Ray{
public:
    Tuples origin;
    Tuples direction;

    Ray(Tuples _origin, Tuples _direction):origin(_origin),direction(_direction){};

    Tuples position(float t)const{
        return origin+direction*t;
    }
};

class Intersection{
public:
    float t;
    const Sphere* object;
    Intersection(float _t, const Sphere* sph): t(_t), object(sph){};
};

class Intersections{
public:
    Intersections(initializer_list<Intersection> intersections): intersects(intersections){};

    Intersection operator[](int index) const{
        return intersects[index];
    }

    auto begin() { return intersects.begin(); }
    auto end()   { return intersects.end(); }

    auto begin() const { return intersects.cbegin(); }
    auto end()   const { return intersects.cend(); }

    int count()const{
        return intersects.size();
    }
private:
    vector<Intersection> intersects;
};

Ray transform(Ray& ab, fourmatrix& mat){
    Tuples orig = mat*ab.origin;
    Tuples dir = mat*ab.direction;
    return Ray(orig, dir);
}

Intersections Intersects(Ray& ray,Sphere& spheres){
    fourmatrix tran= spheres.tranform.invertible();
    Ray rays = transform(ray,tran);
    Tuples sphere_to_ray = rays.origin - point(0,0,0);
    float a = rays.direction * rays.direction;
    float b = 2* (rays.direction * sphere_to_ray);
    float c = (sphere_to_ray * sphere_to_ray) -1;
    float discriminant = pow(b,2)-4*a*c;
    if(discriminant<0){
        return {};
    }
    float t1=(-b-sqrt(discriminant))/(2*a);
    float t2=(-b+sqrt(discriminant))/(2*a);
    return {Intersection(t1,&spheres),Intersection(t2,&spheres)};
};

Intersection hit(Intersections& inters){
    vector<float> ts;
    for (Intersection& a: inters){
        ts.push_back(a.t);
    }
    float smallest=numeric_limits<float>::max();
    int count = 0;
    int index = -1;
    for (float small: ts){
        if(small>0 && small<smallest){
            smallest = small;
            index = count;
        }
        count++;
    }
    if (index!=-1){
        return inters[index];
    }
    return Intersection(0,nullptr);
}

void twelve_hour_clock(){
    const double PI = std::acos(-1.0);
    
    Canvas canva(500,500);
    Pixels white(1,1,1);
    Tuples first = point(0,1,0)*200;

    canva.write_pixels(first.x+250,first.y+250,white);

    Tuples* ptr;
    ptr = &first;
    fourmatrix thirty = rotation_z(PI/6);
    for(int x=0;x<11;x++){
        Tuples new_point = thirty * *ptr;
        cout<<new_point.x<<" "<<new_point.y<<" ";
        canva.write_pixels(new_point.x+250,new_point.y+250,white);
        *ptr = new_point;
    }

    string out = canva.ppm_data();
    // cout<<out;
    ofstream outfile("clock.ppm");
    if (outfile.is_open()){
        outfile<<out;
        outfile.close();
        cout<<"Closed successfully"<<endl;
    }
    else{
        cout<<"Couldn't open successfully";
    }
}

void red_sphere(){
    int canvas_pixels = 100;
    float wall_size = 7;
    float pixel_size = wall_size/canvas_pixels;
    float wall_z = 10;
    Tuples ray_origin = point(0,0,-5);
    float half = wall_size/2;
    Canvas canva(canvas_pixels,canvas_pixels);
    Pixels red(1,0,0);
    Sphere s(1);
    s.tranform.print();
    for (int y=0;y<canvas_pixels;y++){
        float y_position = half -pixel_size*y;
        for(int x = 0;x<canvas_pixels;x++){
            float x_position = -half+pixel_size*x;
            Tuples position = point(x_position,y_position,wall_z);
            Tuples ray_direction = (position - ray_origin);
            Ray r(ray_origin,ray_direction.normalize());
            Intersections xs = Intersects(r,s);
            Intersection hit_check = hit(xs);
            if (hit_check.object)
            {
                canva.write_pixels(x,y,red);
            }
        }
    }
    string out = canva.ppm_data();
    ofstream outfile("sphere.ppm");
    if (outfile.is_open()){
        outfile<<out;
        outfile.close();
        cout<<"Closed successfully"<<endl;
    }
    else{
        cout<<"Couldn't open successfully";
    }
}

int main() {
    cout<<"Try calling the twelve_hour_clock and red_sphere functions"<<endl;
}
