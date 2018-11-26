#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"

inline bool isWellReconstructed(recob::Track track, simb::MCParticle mcp){

  float tsy = track.Start().Y();
  float tey = track.End().Y();
  float tsz = track.Start().Z();
  float tez = track.End().Z();

  float msy = mcp.Vy();
  float mey = mcp.EndY();
  float msz = mcp.Vz();
  float mez = mcp.EndZ();

  float twoDStartRes = std::sqrt(std::pow(msy-tsy,2)+std::pow(msz-tsz,2));
  float twoDStartResFlip = std::sqrt(std::pow(msy-tey,2)+std::pow(msz-tez,2));
  float twoDEndRes = std::sqrt(std::pow(mey-tey,2)+std::pow(mez-tez,2));
  float twoDEndResFlip = std::sqrt(std::pow(mey-tsy,2)+std::pow(mez-tsz,2));

  if ((twoDStartRes < 2.0 && twoDEndRes < 2.0) || 
      (twoDStartResFlip < 2.0 && twoDEndResFlip < 2.0)){
    std::cout << "[MCP Matching] Found a match!" << std::endl;
    std::cout << "[MCP Matching] Start Res Fwd : " << twoDStartRes << std::endl;
    std::cout << "[MCP Matching] End Res Fwd   : " << twoDEndRes << std::endl;
    std::cout << "[MCP Matching] Start Res Bwd : " << twoDStartResFlip << std::endl;
    std::cout << "[MCP Matching] End Res Bwd   : " << twoDEndResFlip << std::endl;
    return true;
  }
  else {
    /*    std::cout << "[MCP Matching] Start Res Fwd : " << twoDStartRes << std::endl;
          std::cout << "[MCP Matching] End Res Fwd   : " << twoDEndRes << std::endl;
          std::cout << "[MCP Matching] Start Res Bwd : " << twoDStartResFlip << std::endl;
          std::cout << "[MCP Matching] End Res Bwd   : " << twoDEndResFlip << std::endl;
          */    return false;
  }
}
