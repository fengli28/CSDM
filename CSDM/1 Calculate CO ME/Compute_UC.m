function UC=Compute_UC(A,B)
% calculate the uncertainty coefficent for two binary vectors
% A, B: two binary vectors
% UC: the uncertainty coefficent of (A,B)
% Note that the UC(A,B) is not equal to UC(B,A)

ea=EntropyBinaryVector(A); % the entroy of vector A
eb=EntropyBinaryVector(B); % the entroy of vector B
Ie=JointEntropyBinaryVector(A,B); % the joint entroy of vector A and B
if eb==0
    UC=0;
else
    UC=(ea+eb-Ie)/eb;
end
