// points, vectors, frames in 3D

/***************************************/
/* class Vec3 */
class Vec3 { 

  float x, y, z;
  
  Vec3 () { }
  
  Vec3 (float px, float py, float pz) {
    x = px; 
    y = py; 
    z = pz;
  }
  
  Vec3 (float px, float py) {
    x = px; 
    y = py;
  }

  Vec3 ()
  
  Vec3 set (float px, float py, float pz) {
    x = px; 
    y = py; 
    z = pz; 
    return this;
  }; 
  Vec3 setTo(Vec3 V) {
    x = V.x; 
    y = V.y; 
    z = V.z; 
    return this;
  }; 
  Vec3 set (Vec3 V) {
    x = V.x; 
    y = V.y; 
    z = V.z; 
    return this;
  }; 
  Vec3 add(Vec3 V) {
    x+=V.x; 
    y+=V.y; 
    z+=V.z; 
    return this;
  };
  Vec3 add(float s, Vec3 V) {
    x+=s*V.x; 
    y+=s*V.y; 
    z+=s*V.z; 
    return this;
  };
  Vec3 sub(Vec3 V) {
    x-=V.x; 
    y-=V.y; 
    z-=V.z; 
    return this;
  };
  Vec3 mul(float f) {
    x*=f; 
    y*=f; 
    z*=f; 
    return this;
  };
  Vec3 div(float f) {
    x/=f; 
    y/=f; 
    z/=f; 
    return this;
  };
  Vec3 div(int f) {
    x/=f; 
    y/=f; 
    z/=f; 
    return this;
  };
  Vec3 rev() {
    x=-x; 
    y=-y; 
    z=-z; 
    return this;
  };
  float norm() {
    return(sqrt(sq(x)+sq(y)+sq(z)));
  }; 
  Vec3 normalize() {
    float n=norm(); 
    if (n>0.000001) {
      div(n);
    }; 
    return this;
  };
  Vec3 rotate(float a, Vec3 I, Vec3 J) { // Rotate this by angle a parallel in plane (I,J) Assumes I and J are orthogonal
    float x=d(this, I), y=d(this, J); // dot products
    float c=cos(a), s=sin(a); 
    add(x*c-x-y*s, I); 
    add(x*s+y*c-y, J); 
    return this;
  };
} // end class Vec3

/***************************************/
/* clss pt */
class pt { 
  float x=0, y=0, z=0; 
  pt () {
  }; 
  pt (float px, float py) {
    x = px; 
    y = py;
  };
  pt (float px, float py, float pz) {
    x = px; 
    y = py; 
    z = pz;
  };
  pt set (float px, float py, float pz) {
    x = px; 
    y = py; 
    z = pz; 
    return this;
  }; 
  pt set (pt P) {
    x = P.x; 
    y = P.y; 
    z = P.z; 
    return this;
  }; 
  pt setTo(pt P) {
    x = P.x; 
    y = P.y; 
    z = P.z; 
    return this;
  }; 
  pt setTo(float px, float py, float pz) {
    x = px; 
    y = py; 
    z = pz; 
    return this;
  }; 
  pt add(pt P) {
    x+=P.x; 
    y+=P.y; 
    z+=P.z; 
    return this;
  };
  pt add(Vec3 V) {
    x+=V.x; 
    y+=V.y; 
    z+=V.z; 
    return this;
  };
  pt sub(Vec3 V) {
    x-=V.x; 
    y-=V.y; 
    z-=V.z; 
    return this;
  };
  pt add(float s, Vec3 V) {
    x+=s*V.x; 
    y+=s*V.y; 
    z+=s*V.z; 
    return this;
  };
  pt sub(pt P) {
    x-=P.x; 
    y-=P.y; 
    z-=P.z; 
    return this;
  };
  pt mul(float f) {
    x*=f; 
    y*=f; 
    z*=f; 
    return this;
  };
  pt div(float f) {
    x/=f; 
    y/=f; 
    z/=f; 
    return this;
  };
  pt div(int f) {
    x/=f; 
    y/=f; 
    z/=f; 
    return this;
  };
}

/*************************************************************/
// =====  Vec3tor functions
Vec3 V() {
  return new Vec3();
};                                                                          // make Vec3tor (x,y,z)
Vec3 V(float x, float y, float z) {
  return new Vec3(x, y, z);
};                                            // make Vec3tor (x,y,z)
Vec3 V(Vec3 V) {
  return new Vec3(V.x, V.y, V.z);
};                                                          // make copy of Vec3tor V
Vec3 A(Vec3 A, Vec3 B) {
  return new Vec3(A.x+B.x, A.y+B.y, A.z+B.z);
};                                       // A+B
Vec3 A(Vec3 U, float s, Vec3 V) {
  return V(U.x+s*V.x, U.y+s*V.y, U.z+s*V.z);
};                               // U+sV
Vec3 M(Vec3 U, Vec3 V) {
  return V(U.x-V.x, U.y-V.y, U.z-V.z);
};                                              // U-V
Vec3 M(Vec3 V) {
  return V(-V.x, -V.y, -V.z);
};                                                              // -V
Vec3 V(Vec3 A, Vec3 B) {
  return new Vec3((A.x+B.x)/2.0, (A.y+B.y)/2.0, (A.z+B.z)/2.0);
}                      // (A+B)/2
Vec3 V(Vec3 A, float s, Vec3 B) {
  return new Vec3(A.x+s*(B.x-A.x), A.y+s*(B.y-A.y), A.z+s*(B.z-A.z));
};      // (1-s)A+sB
Vec3 V(Vec3 A, Vec3 B, Vec3 C) {
  return new Vec3((A.x+B.x+C.x)/3.0, (A.y+B.y+C.y)/3.0, (A.z+B.z+C.z)/3.0);
};  // (A+B+C)/3
Vec3 V(Vec3 A, Vec3 B, Vec3 C, Vec3 D) {
  return V(V(A, B), V(C, D));
};                                         // (A+B+C+D)/4
Vec3 V(float s, Vec3 A) {
  return new Vec3(s*A.x, s*A.y, s*A.z);
};                                           // sA
Vec3 V(float a, Vec3 A, float b, Vec3 B) {
  return A(V(a, A), V(b, B));
}                                       // aA+bB 
Vec3 V(float a, Vec3 A, float b, Vec3 B, float c, Vec3 C) {
  return A(V(a, A, b, B), V(c, C));
}                   // aA+bB+cC
Vec3 V(pt P, pt Q) {
  return new Vec3(Q.x-P.x, Q.y-P.y, Q.z-P.z);
};                                          // PQ
Vec3 U(Vec3 V) {
  float n = V.norm(); 
  if (n<0.0000001) return V(0, 0, 0); 
  else return V.div(n);
};             // V/||V||
Vec3 U(pt P, pt Q) {
  return U(V(P, Q));
};                                                                 // PQ/||PQ||
Vec3 U(float x, float y, float z) {
  return U(V(x, y, z));
};                                               // make Vec3tor (x,y,z)
Vec3 N(Vec3 U, Vec3 V) {
  return V( U.y*V.z-U.z*V.y, U.z*V.x-U.x*V.z, U.x*V.y-U.y*V.x);
};                  // UxV cross product (normal to both)
Vec3 N(pt A, pt B, pt C) {
  return N(V(A, B), V(A, C));
};                                                   // normal to triangle (A,B,C), not normalized (proportional to area)
Vec3 B(Vec3 U, Vec3 V) {
  return U(N(N(U, V), U));
}        
Vec3 Normal(Vec3 V) {
  if (abs(V.z)<=min(abs(V.x), abs(V.y))) return V(-V.y, V.x, 0); 
  if (abs(V.x)<=min(abs(V.z), abs(V.y))) return V(0, -V.z, V.y);
  return V(V.z, 0, -V.x);
}


/*************************************************************/
// ===== point functions
pt P() {
  return new pt();
};                                                                          // point (x,y,z)
pt P(float x, float y, float z) {
  return new pt(x, y, z);
};                                            // point (x,y,z)
pt P(float x, float y) {
  return new pt(x, y);
};                                                       // make point (x,y)
pt P(pt A) {
  return new pt(A.x, A.y, A.z);
};                                                           // copy of point P
pt P(pt A, float s, pt B) {
  return new pt(A.x+s*(B.x-A.x), A.y+s*(B.y-A.y), A.z+s*(B.z-A.z));
};        // A+sAB
pt L(pt A, float s, pt B) {
  return new pt(A.x+s*(B.x-A.x), A.y+s*(B.y-A.y), A.z+s*(B.z-A.z));
};        // A+sAB
pt P(pt A, pt B) {
  return P((A.x+B.x)/2.0, (A.y+B.y)/2.0, (A.z+B.z)/2.0);
}                             // (A+B)/2
pt P(pt A, pt B, pt C) {
  return new pt((A.x+B.x+C.x)/3.0, (A.y+B.y+C.y)/3.0, (A.z+B.z+C.z)/3.0);
};     // (A+B+C)/3
pt P(pt A, pt B, pt C, pt D) {
  return P(P(A, B), P(C, D));
};                                            // (A+B+C+D)/4
pt P(float s, pt A) {
  return new pt(s*A.x, s*A.y, s*A.z);
};                                            // sA
pt A(pt A, pt B) {
  return new pt(A.x+B.x, A.y+B.y, A.z+B.z);
};                                         // A+B
pt P(float a, pt A, float b, pt B) {
  return A(P(a, A), P(b, B));
}                                        // aA+bB 
pt P(float a, pt A, float b, pt B, float c, pt C) {
  return A(P(a, A), P(b, B, c, C));
}                     // aA+bB+cC 
pt P(float a, pt A, float b, pt B, float c, pt C, float d, pt D) {
  return A(P(a, A, b, B), P(c, C, d, D));
}   // aA+bB+cC+dD
pt P(pt P, Vec3 V) {
  return new pt(P.x + V.x, P.y + V.y, P.z + V.z);
}                                 // P+V
pt P(pt P, float s, Vec3 V) {
  return new pt(P.x+s*V.x, P.y+s*V.y, P.z+s*V.z);
}                           // P+sV
pt P(pt O, float x, Vec3 I, float y, Vec3 J) {
  return P(O.x+x*I.x+y*J.x, O.y+x*I.y+y*J.y, O.z+x*I.z+y*J.z);
}  // O+xI+yJ
pt P(pt O, float x, Vec3 I, float y, Vec3 J, float z, Vec3 K) {
  return P(O.x+x*I.x+y*J.x+z*K.x, O.y+x*I.y+y*J.y+z*K.y, O.z+x*I.z+y*J.z+z*K.z);
}  // O+xI+yJ+kZ
void makePts(pt[] C) {
  for (int i=0; i<C.length; i++) C[i]=P();
}
pt ToScreen(pt P) {
  return P(screenX(P.x, P.y, P.z), screenY(P.x, P.y, P.z), 0);
}  // O+xI+yJ+kZ
pt ToModel(pt P) {
  return P(modelX(P.x, P.y, P.z), modelY(P.x, P.y, P.z), modelZ(P.x, P.y, P.z));
}  // O+xI+yJ+kZ


/*************************************************************/
// ===== mouse
pt Mouse() {
  return P(mouseX, mouseY, 0);
};                                          // current mouse location
pt Pmouse() {
  return P(pmouseX, pmouseY, 0);
};
Vec3 MouseDrag() {
  return V(mouseX-pmouseX, mouseY-pmouseY, 0);
};                     // Vec3tor representing recent mouse displacement
pt ScreenCenter() {
  return P(width/2, height/2);
}                                                        //  point in center of  canvas


/*************************************************************/
// ===== measures
float d(Vec3 U, Vec3 V) {
  return U.x*V.x+U.y*V.y+U.z*V.z;
};                                            //U*V dot product
float dot(Vec3 U, Vec3 V) {
  return U.x*V.x+U.y*V.y+U.z*V.z;
};                                            //U*V dot product
float det2(Vec3 U, Vec3 V) {
  return -U.y*V.x+U.x*V.y;
};                                       // U|V det product
float det3(Vec3 U, Vec3 V) {
  return sqrt(d(U, U)*d(V, V) - sq(d(U, V)));
};                              // U|V det product (absolute value)
float m(Vec3 U, Vec3 V, Vec3 W) {
  return d(U, N(V, W));
};                              // (UxV)*W  mixed product, determinant
float m(pt E, pt A, pt B, pt C) {
  return m(V(E, A), V(E, B), V(E, C));
}                                    // det (EA EB EC) is >0 when E sees (A,B,C) clockwise
float n2(Vec3 V) {
  return sq(V.x)+sq(V.y)+sq(V.z);
};                                                   // V*V    norm squared
float n(Vec3 V) {
  return sqrt(n2(V));
};                                                                // ||V||  norm
float d(pt P, pt Q) {
  return sqrt(sq(Q.x-P.x)+sq(Q.y-P.y)+sq(Q.z-P.z));
};                            // ||AB|| distance
float area(pt A, pt B, pt C) {
  return n(N(A, B, C))/2;
};                                               // area of triangle 
float volume(pt A, pt B, pt C, pt D) {
  return m(V(A, B), V(A, C), V(A, D))/6;
};                           // volume of tet 
boolean parallel (Vec3 U, Vec3 V) {
  return n(N(U, V))<n(U)*n(V)*0.00001;
}                              // true if U and V are almost parallel
float angle(Vec3 U, Vec3 V) {
  return acos(d(U, V)/n(V)/n(U));
};                                       // angle(U,V)
boolean cw(Vec3 U, Vec3 V, Vec3 W) {
  return m(U, V, W)>0;
};                                               // (UxV)*W>0  U,V,W are clockwise
boolean cw(pt A, pt B, pt C, pt D) {
  return volume(A, B, C, D)>0;
};                                     // tet is oriented so that A sees B, C, D clockwise 
boolean projectsBetween(pt P, pt A, pt B) {
  return dot(V(A, P), V(A, B))>0 && dot(V(B, P), V(B, A))>0 ;
};
float disToLine(pt P, pt A, pt B) {
  return det3(U(A, B), V(A, P));
};
pt projectionOnLine(pt P, pt A, pt B) {
  return P(A, dot(V(A, B), V(A, P))/dot(V(A, B), V(A, B)), V(A, B));
}

/*************************************************************/
// ===== rotate 
Vec3 R(Vec3 V) {
  return V(-V.y, V.x, V.z);
} // rotated 90 degrees in XY plane
pt R(pt P, float a, Vec3 I, Vec3 J, pt G) {
  float x=d(V(G, P), I), y=d(V(G, P), J); 
  float c=cos(a), s=sin(a); 
  return P(P, x*c-x-y*s, I, x*s+y*c-y, J);
}; // Rotated P by a around G in plane (I,J)
Vec3 R(Vec3 V, float a, Vec3 I, Vec3 J) {
  float x=d(V, I), y=d(V, J); 
  float c=cos(a), s=sin(a); 
  return A(V, V(x*c-x-y*s, I, x*s+y*c-y, J));
}; // Rotated V by a parallel to plane (I,J)
pt R(pt Q, pt C, pt P, pt R) { // returns rotated version of Q by angle(CP,CR) parallel to plane (C,P,R)
  Vec3 I0=U(C, P), I1=U(C, R), V=V(C, Q); 
  float c=d(I0, I1), s=sqrt(1.-sq(c)); 
  if (abs(s)<0.00001) return Q;
  Vec3 J0=V(1./s, I1, -c/s, I0);  
  Vec3 J1=V(-s, I0, c, J0);  
  float x=d(V, I0), y=d(V, J0);
  //  stroke(red); show(C,400,I0); stroke(blue); show(C,400,I1); stroke(orange); show(C,400,J0); stroke(magenta); show(C,400,J1); noStroke();
  return P(Q, x, M(I1, I0), y, M(J1, J0));
} 
pt R(pt Q, float a) {
  float dx=Q.x, dy=Q.y, c=cos(a), s=sin(a); 
  return P(c*dx+s*dy, -s*dx+c*dy, Q.z);
};  // Q rotated by angle a around the origin
pt R(pt Q, float a, pt C) {
  float dx=Q.x-C.x, dy=Q.y-C.y, c=cos(a), s=sin(a); 
  return P(C.x+c*dx-s*dy, C.y+s*dx+c*dy, Q.z);
};  // Q rotated by angle a around point P


/*************************************************************/
// ===== render
void normal(Vec3 V) {
  normal(V.x, V.y, V.z);
};                                          // changes normal for smooth shading
void vertex(pt P) {
  vertex(P.x, P.y, P.z);
};                                           // vertex for shading or drawing
void v(pt P) {
  vertex(P.x, P.y, P.z);
};                                           // vertex for shading or drawing
void nv(Vec3 N) {
  normal(N.x, N.y, N.z);
};                                           // vertex for shading or drawing
void vTextured(pt P, float u, float v) {
  vertex(P.x, P.y, P.z, u, v);
};                          // vertex with texture coordinates
void show(pt P, pt Q) {
  line(Q.x, Q.y, Q.z, P.x, P.y, P.z);
};                       // draws edge (P,Q)
void show(pt P, Vec3 V) {
  line(P.x, P.y, P.z, P.x+V.x, P.y+V.y, P.z+V.z);
};          // shows edge from P to P+V
void show(pt P, float d, Vec3 V) {
  line(P.x, P.y, P.z, P.x+d*V.x, P.y+d*V.y, P.z+d*V.z);
}; // shows edge from P to P+dV
void show(pt A, pt B, pt C) {
  beginShape(); 
  vertex(A);
  vertex(B); 
  vertex(C); 
  endShape(CLOSE);
};                      // volume of tet 
void show(pt A, pt B, pt C, pt D) {
  beginShape(); 
  vertex(A); 
  vertex(B); 
  vertex(C); 
  vertex(D); 
  endShape(CLOSE);
};                      // volume of tet 
void show(pt P, float r) {
  pushMatrix(); 
  translate(P.x, P.y, P.z); 
  sphere(r); 
  popMatrix();
}; // render sphere of radius r and center P
void show(pt P, float s, Vec3 I, Vec3 J, Vec3 K) {
  noStroke(); 
  fill(yellow); 
  show(P, 5); 
  stroke(red); 
  show(P, s, I); 
  stroke(green); 
  show(P, s, J); 
  stroke(blue); 
  show(P, s, K);
}; // render sphere of radius r and center P
void show(pt P, String s) {
  text(s, P.x, P.y, P.z);
}; // prints string s in 3D at P
void show(pt P, String s, Vec3 D) {
  text(s, P.x+D.x, P.y+D.y, P.z+D.z);
}; // prints string s in 3D at P+D
void showShadow(pt P, float r) {
  pushMatrix(); 
  translate(P.x, P.y, 0); 
  scale(1, 1, 0.01); 
  sphere(r); 
  popMatrix();
}
void showShadow(pt P, float r, float h) {
  pushMatrix(); 
  translate(P.x, P.y, h); 
  scale(1, 1, 0.01); 
  sphere(r); 
  popMatrix();
}

String toText(Vec3 V) { 
  return "("+nf(V.x, 1, 5)+","+nf(V.y, 1, 5)+","+nf(V.z, 1, 5)+")";
}


/*************************************************************/
// ==== curve
void bezier(pt A, pt B, pt C, pt D) {
  bezier(A.x, A.y, A.z, B.x, B.y, B.z, C.x, C.y, C.z, D.x, D.y, D.z);
} // draws a cubic Bezier curve with control points A, B, C, D
void bezier(pt [] C) {
  bezier(C[0], C[1], C[2], C[3]);
} // draws a cubic Bezier curve with control points A, B, C, D
pt bezierPoint(pt[] C, float t) {
  return P(bezierPoint(C[0].x, C[1].x, C[2].x, C[3].x, t), bezierPoint(C[0].y, C[1].y, C[2].y, C[3].y, t), bezierPoint(C[0].z, C[1].z, C[2].z, C[3].z, t));
}
Vec3 bezierTangent(pt[] C, float t) {
  return V(bezierTangent(C[0].x, C[1].x, C[2].x, C[3].x, t), bezierTangent(C[0].y, C[1].y, C[2].y, C[3].y, t), bezierTangent(C[0].z, C[1].z, C[2].z, C[3].z, t));
}
void PT(pt P0, Vec3 T0, pt P1, Vec3 T1) {
  float d=d(P0, P1)/3;  
  bezier(P0, P(P0, -d, U(T0)), P(P1, -d, U(T1)), P1);
} // draws cubic Bezier interpolating  (P0,T0) and  (P1,T1) 
void PTtoBezier(pt P0, Vec3 T0, pt P1, Vec3 T1, pt [] C) {
  float d=d(P0, P1)/3;  
  C[0].set(P0); 
  C[1].set(P(P0, -d, U(T0))); 
  C[2].set(P(P1, -d, U(T1))); 
  C[3].set(P1);
} // draws cubic Bezier interpolating  (P0,T0) and  (P1,T1) 
Vec3 Vec3ToCubic (pt A, pt B, pt C, pt D, pt E) {
  return V( (-A.x+4*B.x-6*C.x+4*D.x-E.x)/6, (-A.y+4*B.y-6*C.y+4*D.y-E.y)/6, (-A.z+4*B.z-6*C.z+4*D.z-E.z)/6);
}
Vec3 Vec3ToProp (pt B, pt C, pt D) {
  float cb=d(C, B);  
  float cd=d(C, D); 
  return V(C, P(B, cb/(cb+cd), D));
};  


/*************************************************************/
// ==== perspective
pt Pers(pt P, float d) { 
  return P(d*P.x/(d+P.z), d*P.y/(d+P.z), d*P.z/(d+P.z) );
};
pt InverserPers(pt P, float d) { 
  return P(d*P.x/(d-P.z), d*P.y/(d-P.z), d*P.z/(d-P.z) );
};


/*************************************************************/
// ==== intersection
boolean intersect(pt P, pt Q, pt A, pt B, pt C, pt X) {
  return intersect(P, V(P, Q), A, B, C, X);
} // if (P,Q) intersects (A,B,C), return true and set X to the intersection point
boolean intersect(pt E, Vec3 T, pt A, pt B, pt C, pt X) { // if ray from E along T intersects triangle (A,B,C), return true and set X to the intersection point
  Vec3 EA=V(E, A), EB=V(E, B), EC=V(E, C), AB=V(A, B), AC=V(A, C); 
  boolean s=cw(EA, EB, EC), sA=cw(T, EB, EC), sB=cw(EA, T, EC), sC=cw(EA, EB, T); 
  if ( (s==sA) && (s==sB) && (s==sC) ) return false;
  float t = m(EA, AC, AB) / m(T, AC, AB);
  X.set(P(E, t, T));
  return true;
}
boolean rayIntersectsTriangle(pt E, Vec3 T, pt A, pt B, pt C) { // true if ray from E with direction T hits triangle (A,B,C)
  Vec3 EA=V(E, A), EB=V(E, B), EC=V(E, C); 
  boolean s=cw(EA, EB, EC), sA=cw(T, EB, EC), sB=cw(EA, T, EC), sC=cw(EA, EB, T); 
  return  (s==sA) && (s==sB) && (s==sC) ;
};
boolean edgeIntersectsTriangle(pt P, pt Q, pt A, pt B, pt C) {
  Vec3 PA=V(P, A), PQ=V(P, Q), PB=V(P, B), PC=V(P, C), QA=V(Q, A), QB=V(Q, B), QC=V(Q, C); 
  boolean p=cw(PA, PB, PC), q=cw(QA, QB, QC), a=cw(PQ, PB, PC), b=cw(PA, PQ, PC), c=cw(PQ, PB, PQ); 
  return (p!=q) && (p==a) && (p==b) && (p==c);
}
float rayParameterToIntersection(pt E, Vec3 T, pt A, pt B, pt C) {
  Vec3 AE=V(A, E), AB=V(A, B), AC=V(A, C); 
  return - m(AE, AC, AB) / m(T, AC, AB);
}

float angleDraggedAround(pt G) {  // returns angle in 2D dragged by the mouse around the screen projection of G
  pt S=P(screenX(G.x, G.y, G.z), screenY(G.x, G.y, G.z), 0);
  Vec3 T=V(S, Pmouse()); 
  Vec3 U=V(S, Mouse());
  return atan2(d(R(U), T), d(U, T));
}

float scaleDraggedFrom(pt G) {
  pt S=P(screenX(G.x, G.y, G.z), screenY(G.x, G.y, G.z), 0); 
  return d(S, Mouse())/d(S, Pmouse());
}

// FANS, CONES, AND ARROWS
void disk(pt P, Vec3 V, float r) {  
  Vec3 I = U(Normal(V));
  Vec3 J = U(N(I, V));
  disk(P, I, J, r);
}

void disk(pt P, Vec3 I, Vec3 J, float r) {
  float da = TWO_PI/36;
  beginShape(TRIANGLE_FAN);
  v(P);
  for (float a=0; a<=TWO_PI+da; a+=da) v(P(P, r*cos(a), I, r*sin(a), J));
  endShape();
}


void fan(pt P, Vec3 V, float r) {  
  Vec3 I = U(Normal(V));
  Vec3 J = U(N(I, V));
  fan(P, V, I, J, r);
}

void fan(pt P, Vec3 V, Vec3 I, Vec3 J, float r) {
  float da = TWO_PI/36;
  beginShape(TRIANGLE_FAN);
  v(P(P, V));
  for (float a=0; a<=TWO_PI+da; a+=da) v(P(P, r*cos(a), I, r*sin(a), J));
  endShape();
}

void collar(pt P, Vec3 V, float r, float rd) {
  Vec3 I = U(Normal(V));
  Vec3 J = U(N(I, V));
  collar(P, V, I, J, r, rd);
}

void collar(pt P, Vec3 V, Vec3 I, Vec3 J, float r, float rd) {
  float da = TWO_PI/36;
  beginShape(QUAD_STRIP);
  for (float a=0; a<=TWO_PI+da; a+=da) {
    v(P(P, r*cos(a), I, r*sin(a), J, 0, V)); 
    v(P(P, rd*cos(a), I, rd*sin(a), J, 1, V));
  }
  endShape();
}

void cone(pt P, Vec3 V, float r) {
  fan(P, V, r); 
  disk(P, V, r);
}

void stub(pt P, Vec3 V, float r, float rd) {
  collar(P, V, r, rd); 
  disk(P, V, r); 
  disk(P(P, V), V, rd);
}

void arrow(pt P, Vec3 V, float r) {
  stub(P, V(.8, V), r*2/3, r/3); 
  cone(P(P, V(.8, V)), V(.2, V), r);
}  

// **************************** PRIMITIVE
void showFrame(float d) { 
  noStroke(); 
  fill(metal); 
  sphere(d/10);
  fill(blue);  
  showArrow(d, d/10);
  fill(red); 
  pushMatrix(); 
  rotateY(PI/2); 
  showArrow(d, d/10); 
  popMatrix();
  fill(green); 
  pushMatrix(); 
  rotateX(-PI/2); 
  showArrow(d, d/10); 
  popMatrix();
}

void showFan(float d, float r) {
  float da = TWO_PI/36;
  beginShape(TRIANGLE_FAN);
  vertex(0, 0, d);
  for (float a=0; a<=TWO_PI+da; a+=da) vertex(r*cos(a), r*sin(a), 0);
  endShape();
}

void showCollar(float d, float r, float rd) {
  float da = TWO_PI/36;
  beginShape(QUAD_STRIP);
  for (float a=0; a<=TWO_PI+da; a+=da) {
    vertex(r*cos(a), r*sin(a), 0); 
    vertex(rd*cos(a), rd*sin(a), d);
  }
  endShape();
}

void showCone(float d, float r) {
  showFan(d, r);  
  showFan(0, r);
}

void showStub(float d, float r, float rd) {
  showCollar(d, r, rd); 
  showFan(0, r);  
  pushMatrix(); 
  translate(0, 0, d); 
  showFan(0, rd); 
  popMatrix();
}

void showArrow(float d, float r) {
  float dd=d/5;
  showStub(d-dd, r*2/3, r/3); 
  pushMatrix(); 
  translate(0, 0, d-dd); 
  showCone(dd, r); 
  popMatrix();
}  

void showBlock(float w, float d, float h, float x, float y, float z, float a) {
  pushMatrix(); 
  translate(x, y, z); 
  rotateZ(TWO_PI*a); 
  box(w, d, h); 
  popMatrix();
}


/*************************************************************/
//*********** PICK
Vec3 I=V(1, 0, 0), J=V(0, 1, 0), K=V(0, 0, 1); // screen projetions of global model frame

void computeProjectedVec3tors() { 
  pt O = ToScreen(P(0, 0, 0));
  pt A = ToScreen(P(1, 0, 0));
  pt B = ToScreen(P(0, 1, 0));
  pt C = ToScreen(P(0, 0, 1));
  I=V(O, A);
  J=V(O, B);
  K=V(O, C);
}

Vec3 ToIJ(Vec3 V) {
  float x = det2(V, J) / det2(I, J);
  float y = det2(V, I) / det2(J, I);
  return V(x, y, 0);
}

Vec3 ToK(Vec3 V) {
  float z = dot(V, K) / dot(K, K);
  return V(0, 0, z);
}
