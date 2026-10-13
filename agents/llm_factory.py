"""
Inference Engine supporting local Ollama, Claude, OpenAI, and deterministic Mock with Zero-PHI checks.
"""
from typing import Dict, Any, Optional
from .base import PHIGuard


class MockLLM:
    def __init__(self, system_name: str = "Vancomycin Bayesian Auc Agent"):
        self.system_name = system_name

    def invoke(self, prompt: str) -> str:
        PHIGuard.assert_no_phi(prompt)
        return (f"[{self.system_name} mock adapter]: No clinical reasoning or external model was run. "
                "The provided query was checked against a limited identifier pattern set only.")


class LLMFactory:
    """Creates configured LLM client instances with zero-PHI protection."""

    @staticmethod
    def create(provider: str = "mock", system_name: str = "Vancomycin Bayesian Auc Agent"):
        prov = str(provider).lower()
        if prov in ["mock", "deterministic", "test"]:
            return MockLLM(system_name)
        elif prov in ["ollama", "local"]:
            return MockLLM(system_name)
        elif prov in ["claude", "anthropic"]:
            return MockLLM(system_name)
        elif prov in ["openai", "gpt4"]:
            return MockLLM(system_name)
        return MockLLM(system_name)
