{
  description = "YOUR DESCRIPTION HERE";

  inputs = {
    # grab nixpkgs, I use unstable!
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";

    # for 'foreach' system
    utils.url = "github:numtide/flake-utils";

    # grab zig overlay for zig
    zig-flake.url = "github:mitchellh/zig-overlay";

    # put our zig into zls to ensure it matches
    zls-flake = {
      url = "github:zigtools/zls?ref=0.15.0";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.zig-overlay.follows = "zig-flake";
    };
  };

  outputs = { self, nixpkgs, utils, zig-flake, zls-flake }:
    utils.lib.eachSystem [ "x86_64-linux" ] (system:
      let

        # packages for the given system
        pkgs = import nixpkgs {
          inherit system;
          # use overlays
          overlays = [
            (final: prev: {
              zig = zig-flake.packages.${system}."0.15.1";
              zls = zls-flake.packages.${system}.default.overrideAttrs (old: {
                nativeBuildInputs = (old.nativeBuildInputs or [ ])
                  ++ [ final.zig ];
              });
              fftwNative =
                (prev.fftw.override { precision = "single"; }).overrideAttrs
                (old: {
                  allowSubstitutes = false;
                  preferLocalBuild = true;
                  NIX_CFLAGS_COMPILE = (old.NIX_CFLAGS_COMPILE or "")
                    + "-march=skylake-avx512 -mtune=skylake-avx512";
                  configureFlags = (old.configureFlags or [ ])
                    ++ (if final.stdenv.isx86_64 then [
                      "--enable-sse2"
                      "--enable-avx"
                      "--enable-avx2"
                      "--enable-avx512"
                    ] else
                      [ ]);
                });
            })
          ];
        };
      in {
        # on `nix develop`
        devShells.default = pkgs.mkShell {
          nativeBuildInputs = [ pkgs.zig pkgs.zls pkgs.fftwNative ];

          # puts a nice hook, I like this
          shellHook = ''
            PS1="(dev) $PS1"
          '';
        };
      });
}
