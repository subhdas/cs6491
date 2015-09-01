class Collison{

  int N;
  float[] Ball2Ball, Ball2Wall;

  Collison(int _N, Balls balls){
    N = _N;
    Ball2Ball = new float[N * (N-1) / 2];
    Ball2Wall = new float[N * 6];

    for( int i = 0; i < N - 1; i++ ){
      for(int j = i+1; j < N; j++){
        
      }
    }
  }

  // calculate the collison time between two balls
  // given their locations and velocities.
  void calCollisionTime(Vec3 pa, Vec3 va, Vec3 pb, Vec3 vb){
    Vec3 ab = Vec3(b.sub(a));
    va
  }

}
