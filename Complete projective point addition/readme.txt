%----------------------------------------------------------------------
%-- Author: Poncha Joel
%--
%-- Project Name: Point Addition in Eliptic curve cryptography
%-- Module Name: projective point addition.
%-- Description: This is a matlab implementation of the Compleate projective point addition for prime order a=0
% short Weistrass curve y^2 = x^3 + b. The equations are proposed in this
% work "Complete addition formulas for prime order elliptic curves"  https://eprint.iacr.org/2015/1060.pdf
%----------------------------------------------------------------------
% Remark: Use this website to select curve and points and verify the
% addition results in this matlab program
% (https://andrea.corbellini.name/ecc/interactive/modk-add.html).
% This site can validate addition for big numbers: http://www.christelbach.com/eccalculator.aspx
%----------------------------------------------------------------------
             HOW TO EXECUTE
%-----------------------------------------------------------------------
  1. Download this repo with all files.
  2. Open folder in MATLAB
  3. Run root.m and modify it according to directions in comments
  
  Note: for SECP256K1, the ouput shown in the comments (X3,Y3,Z3) is what is expected, but the equations do not achieve this.
