#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/massmatrix.h>
#include <igl/readOFF.h>
#include <iostream>

#define TUTORIAL_SHARED_PATH "/home/davidbelgrod/libigl/tutorial/data"

Eigen::MatrixXd V; // =(Eigen::MatrixXd(8,3)<<
//     0.0,0.0,0.0,
//     0.0,0.0,1.0,
//     0.0,1.0,0.0,
//     0.0,1.0,1.0,
//     1.0,0.0,0.0,
//     1.0,0.0,1.0,
//     1.0,1.0,0.0,
//     1.0,1.0,1.0).finished();
Eigen::MatrixXi F, T;// = (Eigen::MatrixXi(12,3)<<
    // 1,7,5,
    // 1,3,7,
    // 1,4,3,
    // 1,2,4,
    // 3,8,7,
    // 3,4,8,
    // 5,7,8,
    // 5,8,6,
    // 1,5,6,
    // 1,6,2,
    // 2,6,8,
    // 2,8,4).finished().array()-1;

Eigen::MatrixXd U;
//Eigen::MatrixXi F;
Eigen::SparseMatrix<double> L;
igl::opengl::glfw::Viewer viewer;
int iter;
float lambda;


int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  if (argc < 2) {
        cout << "Usage ex1_bin mesh.obj" << endl;
        exit(0);
    }
  
  const std::string word = argv[1];
  iter = 1000;
  lambda = 1/1000;

  if (argc >= 3) {
      iter = atoi(argv[2]);
  }
  
  if (argc >= 4) {
      lambda = atof(argv[3]);
  }
  
  const std::string join = "/";
  // // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH + join +  word, V, F);
 
  // Compute Laplace-Beltrami operator: #V by #V
  igl::cotmatrix(V,F,L);

  // Alternative construction of same Laplacian
  SparseMatrix<double> G,K;
  // Gradient/Divergence
  igl::grad(V,F,G);
  // Diagonal per-triangle "mass matrix"
  VectorXd dblA;
  igl::doublearea(V,F,dblA);
  // Place areas along diagonal #dim times
  const auto & T = 1.*(dblA.replicate(3,1)*0.5).asDiagonal();
  // Laplacian K built as discrete divergence of gradient or equivalently
  // discrete Dirichelet energy Hessian
  K = -G.transpose() * T * G;
  cout<<"|K-L|: "<<(K-L).norm()<<endl;

  const auto &key_down = [](igl::opengl::glfw::Viewer &viewer,unsigned char key,int mod)->bool
  {
    switch(key)
    {
      case 'r':
      case 'R':
        U = V;
        break;
      case ' ':
      {
      #pragma omp parallel for //useless bc we cant parallelize
        for (int j = 0; j < iter; j++){
          cout << "Smoothing iteration: " << j << endl;
          // Recompute just mass matrix on each step
          SparseMatrix<double> M;
          igl::massmatrix(U,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
          
          // implicit smoothing
          // // Solve (M-delta*L) U = M*U
          // const auto & S = (M- lambda*L);
          // Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
          // assert(solver.info() == Eigen::Success);
          // U = solver.solve(M*U).eval();
          
          //explicit smoothing
          const auto & S = (M + lambda *L);
          SimplicialLLT<SparseMatrix<double > > solver(M);
          assert(solver.info() == Success);
          U = solver.solve(S*U).eval();

          // Compute centroid and subtract (also important for numerics)
          VectorXd dblA;
          igl::doublearea(U,F,dblA);
          double area = 0.5*dblA.sum();
          MatrixXd BC;
          igl::barycenter(U,F,BC);
          RowVector3d centroid(0,0,0);
          for(int i = 0;i<BC.rows();i++)
          {
            centroid += 0.5*dblA(i)/area*BC.row(i);
          }
          U.rowwise() -= centroid;
          // Normalize to unit surface area (important for numerics)
          U.array() /= sqrt(area);
          
          float ratio;
          int width, height;
        
 
          // glfwGetFramebufferSize(viewer.window, &width, &height); 
          // //gets width, height and stores them in addresses to fetch later.
          
          // //int* width = means this a pointer 
          // ratio = width / (float) height;
  
          // glViewport(0, 0, width, height);

          //if ((j % 5) == 0) {
            viewer.data().set_vertices(U);
            viewer.data().compute_normals();
            viewer.core.align_camera_center(U,F);
            viewer.draw();
            glfwSwapBuffers(viewer.window);
            glfwPollEvents();
         // }
        //     if (j % 10 == 0){
        //       lambda *= 5;
        //     }
        }
        break;
        
      }
      default:
        return false;
    }
    // Send new positions, update normals, recenter
    viewer.data().set_vertices(U);
    viewer.data().compute_normals();
    viewer.core.align_camera_center(U,F);
    

    return true;
  };

  // Use original normals as pseudo-colors
  MatrixXd N;
  igl::per_vertex_normals(V,F,N);
  MatrixXd C = N.rowwise().normalized().array()*0.5+0.5;

  // Initialize smoothing with base mesh
  U = V;
  viewer.data().set_mesh(U, F);
  viewer.data().set_colors(C); 
  viewer.callback_key_down = key_down;
  viewer.core.is_animating = true;

  cout<<"Press [space] to smooth."<<endl;;
  cout<<"Press [r] to reset."<<endl;;
  return viewer.launch();
}
