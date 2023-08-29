let
  mozillaOverlay =
    import (builtins.fetchGit {
      url = "https://github.com/mozilla/nixpkgs-mozilla.git";
      rev = "78e723925daf5c9e8d0a1837ec27059e61649cb6";
    });
  nixpkgs = import <nixpkgs> { overlays = [ mozillaOverlay ]; };
  rust-nightly = with nixpkgs; ((rustChannelOf { date = "2022-11-15"; channel = "nightly"; }).rust);
in
with nixpkgs; pkgs.mkShell {
  nativeBuildInputs = [
    rust-nightly
  ];
}
