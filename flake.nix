{
  description = "C/C++ FHS development environment for x86_64-linux";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs =
    { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs { inherit system; };

      # Build the FHS environment
      fhs = pkgs.buildFHSUserEnv {
        name = "clang-fhs-env";

        # Packages accessible in standard paths (/usr/bin, /usr/include, /lib)
        targetPkgs =
          pkgs: with pkgs; [
            clang
            llvmPackages.bintools
            glibc.dev
            gnumake
            cmake
            pkg-config

            # Add extra dependencies here
          ];

        runScript = "bash";
      };
    in
    {
      devShells.${system}.default = fhs.env;
    };
}
