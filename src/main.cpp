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
        return Tuples(
            roundTo2Decimals(x/mag),
            roundTo2Decimals(y/mag),
            roundTo2Decimals(z/mag),
            w);
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

    Tuples operator*(const Tuples& ab) const{
        return Tuples(x*ab.x,y*ab.y,z*ab.z,w*ab.w);
    }

    Tuples operator*(const float num) const{
        return Tuples(x*num,y*num,z*num,w*num);
    }

    Tuples operatorX(const Tuples& ab){
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

    Pixels(float _red, float _green, float _blue): 
    red(_red), green(_green), blue(_blue){};

    Pixels operator+(Pixels& c2){
        return Pixels(red+c2.red,green+c2.green, blue+c2.blue);
    }

    Pixels operator-(Pixels& c2){
        return Pixels(red-c2.red,green-c2.green, blue-c2.blue);
    }

    Pixels operator*(const int x){
        return Pixels(red*x,green*x, blue*x);
    }

    Pixels operator*(Pixels& c2){
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

    twomatrix()= default;

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
    threematrix() = default;

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
    fourmatrix() = default;
    
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

    fourmatrix operator*(fourmatrix& b){
        fourmatrix news;
        for(int row=0;row<4;row++){
            for(int col=0;col<4;col++){
                news[row][col]= data[row][0]*b[0][col]+ data[row][1]*b[1][col] + data[row][2]*b[2][col] + data[row][3]*b[3][col];
            }
        }
        return news;
    }

    fourmatrix operator/(float a)const{
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

    Tuples operator*(Tuples& ab){
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

int main() {

    const double PI = std::acos(-1.0);

    // Tuples first = point(2,3,4);
    // fourmatrix ab = shearing(0,0,1,0,0,0);
    // Tuples tranformed = ab * first;
    // tranformed.display();
    
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
    ofstream outfile("image.ppm");
    if (outfile.is_open()){
        outfile<<out;
        outfile.close();
        cout<<"Closed successfully"<<endl;
    }
    else{
        cout<<"Couldn't open successfully";
    }

    // Canvas canva(900,550);
    // const Pixels c1(1.5,0,0);

    // const Tuples grav = vector_dec(0,-0.1,0);
    // const Tuples wind = vector_dec(-0.01,0,0);
    // grav.point_or_vec();
    // wind.point_or_vec();
    // const environment env(grav,wind);

    // Tuples firstpoint = point(0,1,0);
    // Tuples firstvelocity = vector_dec(1,1.8,0).normalize() * 11.25;
    // projectile proj = projectile(firstpoint, firstvelocity);

    // while (proj.position.y>=0){
    //     proj = tick(proj, env);
    //     float x = round(proj.position.x);
    //     float y = canva.height - round(proj.position.y);
    //     canva.write_pixels(x,y,c1);
    // }

    

}
