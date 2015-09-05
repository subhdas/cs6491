class BALLS {          // class for manipulaitng and displaying points
  Boolean loop=true;
  int pv =0,                 // index of picked ball
    iv=0,                    // insertion vertex index
    nv = 0;                  // number of vertices currently used in P
  int maxnv = 16000;         // max number of vertices
  pt[] G = new pt [maxnv];   // ball centers
  vec[] V = new vec [maxnv]; // velocities
  float[] r = new float [maxnv]; // radii
  color[] c = new color [maxnv]; // colors
  int[] m = new int [maxnv];     // delay before resetting color
  
  BALLS() {
  }

  BALLS declare() {
    for (int i=0; i<maxnv; i++) {
      G[i]=P(); 
      V[i]=V(); 
      r[i]=1; 
      c[i]=yellow; 
      m[i]=0;
    }
    return this;
  }     // init all point objects

  /**
   * initialize balls on a grid [nb x nb x nb]
   * @param nb       number of BALLS in each dimension
   * @param w        length of cube
   * @param r        radius of the ball
   * @param c        color of the ball
   * @param vRange   initial velocity range
   */
  void initPointsOnGrid(int nb, float w, float r, color c, float vRange) { 
    empty();
    float d = w / (nb+1), dd=d;
    for (int i=0; i<nb; i++) 
      for (int j=0; j<nb; j++) 
        for (int k=0; k<nb; k++)
          addBall(P(d*i-w/2+dd, d*j-w/2+dd, d*k-w/2+dd),
                  V(random(-vRange, vRange), random(-vRange, vRange), random(-vRange, vRange)), r, c);
  }



  //////////////////////////////////////////////////////////////////////
  //                         add by student                           //
  //////////////////////////////////////////////////////////////////////
  class Colli{
    int p;
    float t;

    Colli (int index, float time){
      p = index;
      t = time;
    }
  }

  // static final int Wxp = -1;    // positive x wall
  // static final int Wxn = -2;    // negative x wall
  // static final int Wyp = -3;
  // static final int Wyn = -4;
  // static final int Wzp = -5;
  // static final int Wzn = -6;

  Colli[] colli;                // nv element array
  int minIndex;                 // next collision index
  float minTime;                // next collision time 
  float internalTime;           // internal timer
  /**
   * this alternative initialization function is used to test the one to one collision
   * for debug purpose.
   */
  // void initPointsOnGrid(int nb, float w, float r, color c, float vRange) { 
  //   empty();

  //   int No = 3;
  //   switch (No) {
      
  //   case 1:                     // test bounce wall
  //     addBall(P(-100, 100, 0), V(1, 0, 0), r, c);
  //     addBall(P(100, 100, 0), V(1, 0, 0), r, c);
  //     break;
      
  //   case 2:                     // test ball bounce
  //     addBall(P(-100, 100, 0), V(1, 0, 0), r, c);
  //     addBall(P(100, 100, 0), V(-1, 0, 0), r, c);
  //     break;

  //   case 3: // test a more complicated case
  //     addBall(P(-100, 0, 100), V(1, 0, -1), r, c);
  //     addBall(P(100, 0, -100), V(-1, 0, 1), r, c);

  //   default:
  //     break;
  //   }
  // }


  /**
   * Calculate the collision time between two balls. If two balls cannot
   * collide in positive finite time, the return -1.
   * 
   * @param A    ball A
   * @param B    ball B
   * @param va   the velocity of A
   * @param vb   the velocity of B
   * @param ra   the radius of A
   * @param rb   the radius of B
   */
  float calB2BTime(pt A, pt B, vec va, vec vb, float ra, float rb){
    // solve (A-B)^2 + 2*(A-B)*(VA-VB)*t + (VA-VB)^2*t^2 = (rA+rB)^2
    // => a*t^2 + b*t + c = 0
    vec ab = V(B, A);           // caution ! this is A - B
    vec vab = va.sub2(vb); 

    float a = dot(vab, vab);
    float b = 2 * dot(ab, vab);
    float c = dot(ab, ab) - sq(ra + rb);
    
    if(a == 0) {
      return -1.0;
    }
    else {
      float delta = sq(b) - 4 * a * c;
      if (delta < 0) {
        return -1.0;
      }
      else {
        float t1 = (-b + sqrt(delta)) / (2 * a);
        float t2 = (-b - sqrt(delta)) / (2 * a);
        if (t1 >=0 && t2 >=0) return (t1 > t2 ? t2 : t1);
        else if (t1 >=0 && t2 <=0) return  t1;
        else if (t1 <=0 && t2 >=0) return t2;
        else return -1.0;
      }          
    }        
  }


  /**
   * Calculate the collision time between one ball and one wall.
   * If they cannot collide in positive finite time, then return -1.
   * 
   * @param wall   one dimensional position of wall
   * @param ball   one dimensional position of ball
   * @param ve     one dimensional velocity of ball
   */
  float calB2WTime(float wall, float ball, float v){
    if (v == 0 ) return -1;
    else {
      float time = (wall - ball) / v;
      if (time < 0) return -1;
      else return (wall - ball) / v;
    }
  }

  
  /**
   * Calculate the collision between a ball and the six
   * walls.
   *
   * @param a  the position of the ball
   * @param v  the velocity of the ball
   * @param w  the width of the wall
   * @param r  the radius of the ball
   */
  float[] calB2WTimeAll(pt a, vec v, float w, float r){
    float[] time = new float[6];
    
    time[0] = calB2WTime(-w/2+r, a.x, v.x);
    time[1] = calB2WTime(w/2-r, a.x, v.x);
      
    time[2] = calB2WTime(-w/2+r, a.y, v.y);
    time[3] = calB2WTime(w/2-r, a.y, v.y);

    time[4] = calB2WTime(-w/2+r, a.z, v.z);
    time[5] = calB2WTime(w/2-r, a.z, v.z);

    // for( int i = 0; i < 6; i++) println(time[i]);
    return time;
  }


  /**
   * calculate the smallest collision time between ball i and other
   * balls and walls
   *
   * @param i  the index of one ball
   */
  Colli calOneBallColli(int i) {
    float time = Float.POSITIVE_INFINITY;
    int index = 0;

    float[] t = calB2WTimeAll(G[i], V[i], w, r[i]);
    for (int j = 0; j < 6; j++){
      if(t[j] >= 0 && t[j] < minTime){
        minTime = t[j];
        index = -(j+1);
        }
    }
    
    for (int j = 0; j < nv; j++){
      if (j != i) {
        float t2 = calB2BTime(G[i], G[j], V[i], V[j], r[i], r[j]);
        if (t2 >= 0 && t2 < time) {
          time = t2;
          index = j;
        }
      }
    }
    
    return new Colli(index, time);
  }

  /**
   * Obtain the next collision index ball and update the next collision
   * time.
   */
  int getMinIndex(){
    int minIndex = -1;
    float minT = Float.POSITIVE_INFINITY;
    for (int i = 0; i < colli.length; i++){
      if(colli[i].t < minT){
        minIndex = i;
        minT = colli[i].t;
      }
    }
    return minIndex;
  }
  
  /**
   * Initialize the collision structure.
   * the index of walls are
   *   -1 : negative x wall
   *   -2 : positive x wall
   *   -3 : negative y wall
   *   -4 : positive y wall
   *   -5 : negative z wall
   *   -6 : positive z wall
   */ 
  void initCollision(float w){
    colli = new Colli[nv];
 
    for ( int i = 0; i < nv; i++){
      colli[i] = calOneBallColli(i);
    }

    minIndex = getMinIndex();
    minTime = colli[minIndex].t;
    internalTime = 0;
  }

  void bounceTwoBalls(int i, int j){
    pt A = G[i];
    pt B = G[j];
    vec va = V[i];
    vec vb = V[j];
    
    vec ab = V(A, B);                             // AB
    vec vaParallel = V(dot(va, ab) / n2(ab), ab); // va*ab / |ab|^2 * ab
    vec vaPerpen = M(va, vaParallel);             // perpendicular component
    vec vbParallel = V(dot(vb, ab) / n2(ab), ab);
    vec vbPerpen = M(vb, vbParallel);

    // update the velocity
    V[i] = A(vbParallel, vaPerpen);
    V[j] = A(vaParallel, vbPerpen);
  }

  void bounceBallWall(int i, int wallIndex){
    switch (wallIndex){

    case -1:
    case -2:
      V[i].x *= -1;
      break;
      
    case -3:
    case -4:
      V[i].y *= -1;
      break;

    case -5:
    case -6:
      V[i].z *= -1;
      break;
      
    default:
      println("wrong wall index !");
    }
    
  }

  void updateB2BColli(int k1, int k2) {
    // bounce these two colliding balls
    bounceTwoBalls(k1, k2);
    
    // update the collision time of these two balls
    colli[k1] = calOneBallColli(k1);
    colli[k2] = calOneBallColli(k2);

    // update the collision time of remaining balls
    for( int i = 0; i < nv; i++){
      if (i != k1 && i != k2) {
        // if the minimal time collision is not with k1 or k2
        // this part is linear
        if( colli[i].p != k1 && colli[i].p != k2) {
          float t1 = calB2BTime(G[i], G[k1], V[i], V[k1], r[i], r[k1]);
          float t2 = calB2BTime(G[i], G[k2], V[i], V[k2], r[i], r[k2]);
          float[] ts = {colli[i].t, t1, t2};
          int[] indices = {colli[i].p,  k1, k2};

          float time = Float.POSITIVE_INFINITY;
          int index = 0;
          for (int j = 0; j < 3; j++){
            if(ts[j] < time) {
              time = ts[j];
              index = indices[j];
            }
          }
          colli[i] = new Colli(index, time);
        }
        // if the minimal time collision is with k1 or k2
        else {
          colli[i] = calOneBallColli(i);
        }
      }
    }
    
  }

  void updateB2WColli(int k, int wallIndex){
    // bounce the ball and wall
    bounceBallWall(k, wallIndex);

    // update the collision time of the ball
    colli[k] = calOneBallColli(k);

    // update the collision time of remaining balls
    for( int i = 0; i < nv; i++){
      if (i != k) {
        // if the minimal time collision is not with ball k
        // this part is linear
        if(colli[i].p != k) {
          float t = calB2BTime(G[i], G[k], V[i], V[k], r[i], r[k]);
          if(t < colli[i].t) colli[i] = new Colli(k, t);
        }
        else {
          colli[i] = calOneBallColli(i);
        }
      }
    }
    
  }
  
  void updateColli(){
    int k = colli[minIndex].p;
    
    if (k < 0) updateB2BColli(minIndex, k);
    else updateB2WColli(minIndex, k);

    minIndex = getMinIndex();
    minTime = colli[minIndex].t;
    internalTime = 0;           // reset internal timer
  }

  
  /**
   * @param m  the scaling factor for velocity
   */
  void updateState(float m){
    if (internalTime >= minTime) updateColli();
    internalTime += m;
    advectBalls(m);
  }

  //////////////////////////////////////////////////////////////////////
  //                   END : add by student                           //
  //////////////////////////////////////////////////////////////////////
  
  BALLS addBall(pt Pp, vec Vp, float rp, color cp) { 
    G[nv].setTo(Pp); 
    V[nv].setTo(Vp); 
    r[nv]=rp; 
    c[nv]=cp; 
    pv=nv; 
    nv++;  
    return this;
  } // adds a point at the end

  
  void showBalls() { 
    noStroke(); 
    for (int v=0; v<nv; v++) {
      fill(c[v]); 
      show(G[v], r[v]);         // show sphere at G[v] with radius r[v]
    }
  }

  void showPickedBall() { 
    noStroke(); 
    fill(magenta); 
    show(G[pv], r[pv]+1);
  } 

  void resetColors(color dc) { 
    for (int v=0; v<nv; v++){
      if (m[v]>0) {
        m[v]--; 
        if (m[v]==0) c[v]=dc;
      }
    }
  }
  
  void showVelocities(float m) {
    noStroke(); 
    fill(blue); 
    for (int v=0; v<nv; v++) arrow(G[v], V(m, V[v]), 3);
  }

  void advectBalls(float m) {
    for (int v=0; v<nv; v++) G[v]=P(G[v], V(m, V[v]));
  }

  void bounceBalls(float w) {
    for (int v=0; v<nv; v++) {
      if (G[v].x<-w/2+r[v] || G[v].x>w/2-r[v]) {
        V[v].x=-V[v].x; 
        c[v]=green; 
        m[v]=del;
      }
      if (G[v].y<-w/2+r[v] || G[v].y>w/2-r[v]) {
        V[v].y=-V[v].y; 
        c[v]=green; 
        m[v]=del;
      }
      if (G[v].z<-w/2+r[v] || G[v].z>w/2-r[v]) {
        V[v].z=-V[v].z; 
        c[v]=green; 
        m[v]=del;
      }
    }
    for (int u=0; u<nv-1; u++) for (int v=u+1; v<nv; v++) if (d(G[u], G[v])<r[u]+r[v]) {
      c[u]=red; 
      m[u]=del; 
      c[v]=red; 
      m[v]=del;
    }
  }

  

  BALLS empty() {
    nv=0; 
    pv=0; 
    return this;
  } // resets P so that we can start adding points
  BALLS addPt(pt P) { 
    G[nv].setTo(P); 
    pv=nv; 
    nv++;  
    return this;
  } // adds a point at the end
  BALLS addPt(float x, float y) { 
    G[nv].x=x; 
    G[nv].y=y; 
    pv=nv; 
    nv++; 
    return this;
  }
  BALLS copyFrom(BALLS Q) {
    empty(); 
    nv=Q.nv; 
    for (int v=0; v<nv; v++) G[v]=P(Q.G[v]); 
    return this;
  }
  BALLS setToL(BALLS P, float t, BALLS Q) { // lerp (linear interpolation betwen P and Q
    empty(); 
    nv=min(P.nv, Q.nv); 
    for (int v=0; v<nv; v++) G[v]=L(P.G[v], t, Q.G[v]); 
    return this;
  }
  BALLS resetOnCircle(int k, float r) { // makes new polyloo[p with k  points on a circle around origin
    empty(); // resert P
    pt C = P(); // center of circle
    for (int i=0; i<k; i++) addPt(R(P(C, V(0, -r, 0)), 2.*PI*i/k, C)); // points on z=0 plane
    pv=0; // picked vertex ID is set to 0
    return this;
  } 
  int idOfVertexWithClosestScreenProjectionTo(pt M) { // for picking a vertex with the mouse
    pp=0; 
    for (int i=1; i<nv; i++) if (d(M, ToScreen(G[i]))<=d(M, ToScreen(G[pp]))) pp=i; 
    return pp;
  }

  pt closestProjectionOf(pt M) {   // for picking inserting O. Returns projection but also CHANGES iv !!!!
    pt C = P(G[0]); 
    float d=d(M, C);       
    for (int i=1; i<nv; i++) if (d(M, G[i])<=d) {
      iv=i; 
      C=P(G[i]); 
      d=d(M, C);
    }
    // get adjacent points circularly 
    for (int i=nv-1, j=0; j<nv; i=j++) { 
      pt A = G[i], B = G[j];
      if (projectsBetween(M, A, B) && disToLine(M, A, B)<d) {
        d=disToLine(M, A, B); 
        iv=i; 
        C=projectionOnLine(M, A, B);
      }
    } 
    return C;
  }

  void pickClosestTo(pt M) {   // for picking inserting O. Returns projection but also CHANGES iv !!!!
    pt C = P(G[0]); 
    float d=d(M, C);       
    for (int i=1; i<nv; i++) if (d(M, G[i])<=d) {
      pv=i; 
      C=P(G[i]); 
      d=d(M, C);
    }
  }

  BALLS insertPt(pt P) { // inserts new vertex after vertex with ID iv
    for (int v=nv-1; v>iv; v--) G[v+1].setTo(G[v]); 
    iv++; 
    G[iv].setTo(P);
    nv++; // increments vertex count
    return this;
  }
  BALLS insertClosestProjection(pt M) {  
    pt P = closestProjectionOf(M); // also sets iv
    insertPt(P);
    return this;
  }

  BALLS deletePicked() {
    for (int i=pv; i<nv; i++) G[i].setTo(G[i+1]); 
    pv=max(0, pv-1); 
    nv--;  
    return this;
  }
  BALLS setPt(pt P, int i) { 
    G[i].setTo(P); 
    return this;
  }
  BALLS showPicked() {
    show(G[pv], 13); 
    return this;
  }
  BALLS drawBalls(float r) {
    for (int v=0; v<nv; v++) show(G[v], r); 
    return this;
  }
  BALLS showPicked(float r) {
    show(G[pv], r); 
    return this;
  }
  BALLS drawClosedCurve(float r) {
    for (int v=0; v<nv-1; v++) stub(G[v], V(G[v], G[v+1]), r, r/2);  
    stub(G[nv-1], V(G[nv-1], G[0]), r, r/2); 
    return this;
  }
  BALLS setPickedTo(int pp) {
    pv=pp; 
    return this;
  }
  BALLS movePicked(vec V) { 
    G[pv].add(V); 
    return this;
  }      // moves selected point (index p) by amount mouse moved recently
  BALLS moveAll(vec V) {
    for (int i=0; i<nv; i++) G[i].add(V); 
    return this;
  };   
  pt Picked() {
    return G[pv];
  } 

  void saveBALLS(String fn) {
    String [] inpBALLS = new String [nv+1];
    int s=0;
    inpBALLS[s++]=str(nv);
    for (int i=0; i<nv; i++) {
      inpBALLS[s++]=str(G[i].x)+","+str(G[i].y)+","+str(G[i].z);
    }
    saveStrings(fn, inpBALLS);
  };

  void loadBALLS(String fn) {
    println("loading: "+fn); 
    String [] ss = loadStrings(fn);
    String subBALLS;
    int s=0;   
    int comma, comma1, comma2;   
    float x, y;   
    int a, b, c;
    nv = int(ss[s++]); 
    print("nv="+nv);
    for (int k=0; k<nv; k++) {
      int i=k+s; 
      float [] xy = float(split(ss[i], ",")); 
      G[k].setTo(xy[0], xy[1], xy[2]);
    }
    pv=0;
  };


  
} // end of BALLS class
