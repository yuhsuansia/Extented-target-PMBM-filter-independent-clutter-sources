function meas_in_gate = ellips_gating_clutter(z,stationary_clutter,gating_size)
%ELLIPS_GATING_CLUTTER performs ellipsoidal gating for stationary sources

if isempty(z)
    meas_in_gate = false(0,1);
    return
end

nz = size(z,2);
meas_in_gate = false(nz,1);

Cy = stationary_clutter.X_clutter;

nu = z - repmat(stationary_clutter.x_clutter,[1 nz]);
dist = diag(nu.'*(Cy\nu));

meas_in_gate(dist<gating_size) = true;

end