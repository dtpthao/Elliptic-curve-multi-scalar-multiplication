### Analysis complexity of Joint Signed Digit Representations for ECDSA (Optimize operation: R = kP and R = u1P + u2Q,  P, Q – elliptic curve points):
-	Use Shamir method to decomposite k to n parts (experienced with n from 1 to 5).
-	Convert each part of k to different representation types: multi-dimension binary, joint sparse form (JSF), higher dimensional joint sparse form (d-JSF, d = [2,5]).
-	Compute and minimize the complexity of simultaneous multi-scalar multiplication algorithm using new converted numbers formats by the number of EC Point Double, Additional, and Negative operations. 
-	Optimizes the running time of the algorithms (C/C++, library: miracl). 
